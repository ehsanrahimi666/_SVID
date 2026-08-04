# ---------------------------------------------------------------------------
# Module 7b: physically based screening and synthetic data
# ---------------------------------------------------------------------------

#' Screen indices against radiative-transfer simulations
#'
#' Uses a leaf-canopy radiative transfer model to generate reflectance under
#' controlled variation of one target variable and of confounding variables,
#' then reports how strongly each index responds to the target and how much it
#' is contaminated by the confounders. An index that tracks leaf area index but
#' also tracks soil brightness is not a leaf-area index.
#'
#' This function is an adapter: the physics is supplied by the `ccrtm` or
#' `hsdar` package, which must be installed. Neither is a hard dependency of
#' `svid`. When neither is available, use [svid_simulate()], which generates
#' synthetic data from an explicit empirical model rather than from radiative
#' transfer.
#'
#' @param pool A [svid_pool] or character vector.
#' @param target Name of the target parameter, e.g. `"LAI"` or `"Cab"`.
#' @param confounders Character vector of nuisance parameters to vary.
#' @param ranges Named list of `c(min, max)` ranges for every varied parameter.
#' @param n Number of simulated canopies.
#' @param sensor Sensor whose bands the spectra are resampled to.
#' @param engine `"ccrtm"` or `"hsdar"`.
#' @param seed Random seed.
#'
#' @return A data frame with the target sensitivity, the maximum confounder
#'   sensitivity and their ratio, sorted by the ratio.
#' @export
#' @examples
#' # requires ccrtm or hsdar
#' if (requireNamespace("ccrtm", quietly = TRUE)) {
#'   svid_prosail(c("(N - R) / (N + R)", "(N - S1) / (N + S1)"),
#'                target = "LAI", confounders = "psoil", n = 60)
#' }
svid_prosail <- function(pool, target = "LAI", confounders = c("psoil", "Cab"),
                         ranges = NULL, n = 200L, sensor = "sentinel2",
                         engine = c("ccrtm", "hsdar"), seed = 1L) {
  engine <- match.arg(engine)
  .need(engine)
  if (!inherits(sensor, "svid_sensor")) sensor <- svid_sensor(sensor)
  def <- list(LAI = c(0.1, 7), Cab = c(10, 80), Cw = c(0.005, 0.03),
              Cm = c(0.002, 0.012), psoil = c(0, 1), N = c(1.2, 2.2),
              hspot = c(0.01, 0.5), tts = c(20, 50))
  if (!is.null(ranges)) def <- utils::modifyList(def, ranges)
  vary <- unique(c(target, confounders))
  set.seed(seed)
  pars <- lapply(vary, function(p) stats::runif(n, def[[p]][1], def[[p]][2]))
  names(pars) <- vary
  spec <- .rtm_run(pars, n, engine)
  refl <- .resample_to_sensor(spec$reflectance, spec$wavelength, sensor)
  d <- .svid_data(lapply(seq_len(ncol(refl)), function(j) matrix(refl[, j], n, 1)),
                  pars[[target]],
                  data.frame(id = as.character(seq_len(n)), x = NA, y = NA,
                             group = as.character(seq_len(n)),
                             site = NA, year = NA),
                  "e1", sensor, 1)
  names(d$values) <- colnames(refl); d$bands <- colnames(refl)
  ex <- as.character(pool)
  sens_t <- svid_separability(ex, d, metric = "spearman")
  out <- data.frame(expression = ex,
                    target_sensitivity = sens_t$score[match(ex, sens_t$expression)],
                    stringsAsFactors = FALSE)
  for (cf in setdiff(vary, target)) {
    d2 <- d; d2$response <- pars[[cf]]
    s <- svid_separability(ex, d2, metric = "spearman")
    out[[paste0("conf_", cf)]] <- s$score[match(ex, s$expression)]
  }
  cc <- grep("^conf_", names(out), value = TRUE)
  out$max_confounder <- if (length(cc)) do.call(pmax, out[cc]) else 0
  out$specificity <- out$target_sensitivity / pmax(out$max_confounder, 1e-6)
  out <- out[order(-out$specificity), ]
  rownames(out) <- NULL
  out
}

.rtm_run <- function(pars, n, engine) {
  if (engine == "ccrtm") {
    p <- ccrtm::getDefaults("prosail5")[[1]]$best
    mat <- do.call(rbind, lapply(seq_len(n), function(i) {
      pp <- p
      for (nm in names(pars)) if (nm %in% names(pp)) pp[[nm]] <- pars[[nm]][i]
      as.numeric(ccrtm::fRTM(rho ~ prospect5 + foursail, pars = list(pp))[, 1])
    }))
    list(reflectance = mat, wavelength = 400:2500)
  } else {
    args <- list(TypeLidf = 1)
    for (nm in names(pars)) args[[nm]] <- pars[[nm]]
    sp <- do.call(hsdar::PROSAIL, args)
    list(reflectance = hsdar::spectra(sp), wavelength = hsdar::wavelength(sp))
  }
}

.resample_to_sensor <- function(refl, wl, sensor, fwhm = 40) {
  out <- vapply(sensor$symbols, function(s) {
    c0 <- sensor$wavelength[[s]]
    if (c0 < min(wl) || c0 > max(wl)) return(rep(NA_real_, nrow(refl)))
    w <- exp(-0.5 * ((wl - c0) / (fwhm / 2.355))^2)
    w <- w / sum(w)
    as.numeric(refl %*% w)
  }, numeric(nrow(refl)))
  colnames(out) <- sensor$symbols
  out[, !apply(is.na(out), 2, all), drop = FALSE]
}

#' Synthetic multispectral data for examples and tests
#'
#' Generates reflectance from an explicit empirical mixture model. This is
#' **not** a radiative transfer model and makes no claim to physical accuracy;
#' it exists so that every function in the package has a fast, dependency-free,
#' fully reproducible example, and so that regression workflows can be
#' demonstrated with a known ground truth. Use [svid_prosail()] when physical
#' realism matters.
#'
#' The model is a two-endmember linear mixture. Each sample has a vegetation
#' fraction \eqn{f} and a soil brightness \eqn{s}. Reflectance in band \eqn{b}
#' with central wavelength \eqn{\lambda_b} is
#' \deqn{\rho_b = f \, v(\lambda_b) + (1 - f)\, s \, g(\lambda_b) + \epsilon_b,}
#' where \eqn{v} is a vegetation endmember with chlorophyll absorption in the
#' red, a red-edge rise modelled as a logistic in wavelength, a NIR plateau and
#' water absorption in the SWIR; \eqn{g} is a linearly increasing soil
#' endmember; and \eqn{\epsilon_b} is Gaussian noise. For time series, \eqn{f}
#' follows a double-logistic phenology whose amplitude, start and length differ
#' between classes.
#'
#' @param n Number of samples.
#' @param sensor Sensor name or object supplying the band set.
#' @param epochs Number of acquisition dates.
#' @param response `"class"` for a three-class categorical response, or one of
#'   `"lai"`, `"fcover"`, `"chlorophyll"` for a continuous one.
#' @param classes Class names when `response = "class"`.
#' @param sites Number of simulated sites, written to the `site` column.
#' @param groups Number of parcels.
#' @param noise Standard deviation of the additive noise.
#' @param seed Random seed.
#'
#' @return A [svid_data] object.
#' @export
#' @examples
#' sim <- svid_simulate(n = 100, sensor = "landsat8", epochs = 6)
#' sim
#' reg <- svid_simulate(n = 100, response = "lai")
#' reg
svid_simulate <- function(n = 300L, sensor = "sentinel2", epochs = 8L,
                          response = c("class", "lai", "fcover", "chlorophyll"),
                          classes = c("crop_a", "crop_b", "built"),
                          sites = 1L, groups = 30L, noise = 0.01, seed = 1L) {
  response <- match.arg(response)
  if (!inherits(sensor, "svid_sensor")) sensor <- svid_sensor(sensor)
  set.seed(seed)
  wl <- sensor$wavelength[sensor$symbols]
  veg <- 0.03 + 0.42 / (1 + exp(-(wl - 715) / 18)) -
    0.20 * exp(-0.5 * ((wl - 670) / 30)^2) -
    0.26 / (1 + exp(-(wl - 1500) / 220))
  veg <- pmax(veg + 0.06 * exp(-0.5 * ((wl - 550) / 40)^2), 0.01)
  soil <- 0.10 + 0.00016 * (wl - 400)
  soil <- pmin(pmax(soil, 0.05), 0.55)

  t <- seq(0, 1, length.out = epochs)
  if (response == "class") {
    cls <- factor(sample(classes, n, TRUE), levels = classes)
    par <- data.frame(amp = c(0.85, 0.75, 0.10)[as.integer(cls)],
                      st = c(0.25, 0.45, 0.50)[as.integer(cls)],
                      len = c(0.45, 0.30, 0.90)[as.integer(cls)])
    y <- cls
  } else {
    lvl <- stats::runif(n, 0, 1)
    par <- data.frame(amp = lvl, st = 0.3, len = 0.5)
    y <- switch(response, lai = 7 * lvl, fcover = lvl, chlorophyll = 10 + 70 * lvl)
  }
  bright <- stats::runif(n, 0.6, 1.4)
  vals <- lapply(seq_along(wl), function(b) matrix(NA_real_, n, epochs))
  for (j in seq_len(epochs)) {
    f <- par$amp / (1 + exp(-(t[j] - par$st) / 0.06)) *
      (1 - 1 / (1 + exp(-(t[j] - par$st - par$len) / 0.08)))
    f <- pmin(pmax(f, 0), 1)
    for (b in seq_along(wl)) {
      vals[[b]][, j] <- pmax(f * veg[b] + (1 - f) * bright * soil[b] +
                               stats::rnorm(n, 0, noise), 0.001)
    }
  }
  names(vals) <- sensor$symbols
  meta <- data.frame(id = sprintf("s%04d", seq_len(n)),
                     x = stats::runif(n, 0, 1000), y = stats::runif(n, 0, 1000),
                     group = paste0("g", rep_len(seq_len(groups), n)),
                     site = paste0("site", rep_len(seq_len(sites), n)),
                     year = "2022", stringsAsFactors = FALSE)
  .svid_data(vals, y, meta, sprintf("e%02d", seq_len(epochs)), sensor, 1)
}

#' Small data set used throughout the examples
#'
#' A convenience wrapper around [svid_simulate()] with settings chosen to run
#' in well under a second.
#'
#' @param ... Passed to [svid_simulate()].
#' @return A [svid_data] object.
#' @export
#' @examples
#' svid_demo_data()
svid_demo_data <- function(...) {
  a <- utils::modifyList(list(n = 150L, sensor = "sentinel2", epochs = 5L,
                              groups = 15L, seed = 42L), list(...))
  d <- do.call(svid_simulate, a)
  keep <- intersect(c("B", "G", "R", "N", "S1", "S2"), d$bands)
  d$values <- d$values[keep]; d$bands <- keep
  d
}

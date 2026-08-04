# ---------------------------------------------------------------------------
# Module 7: robustness
# ---------------------------------------------------------------------------

#' Sensitivity of indices to radiometric perturbation
#'
#' A discovered index is only useful if it survives the errors that separate
#' one acquisition from another. This function re-evaluates candidates after
#' perturbing the reflectance data with multiplicative gain, additive offset,
#' band-wise noise and soil-brightness variation, and reports how much each
#' index degrades. Indices whose advantage evaporates under modest perturbation
#' should not be promoted, however well they score on a single scene.
#'
#' @param pool A [svid_pool] or character vector.
#' @param data A [svid_data] object.
#' @param gain,offset,noise,soil Standard deviations of the multiplicative
#'   gain, additive offset (in reflectance units), band-wise Gaussian noise and
#'   soil-brightness perturbation.
#' @param n_rep Number of perturbed replicates.
#' @param metric Screening metric, passed to [svid_separability()].
#' @param seed Random seed.
#'
#' @return A data frame with the unperturbed score, the mean and standard
#'   deviation under perturbation, and the relative loss.
#' @export
#' @examples
#' d <- svid_demo_data()
#' svid_perturb(c("(N - R) / (N + R)", "N / (B * R)", "N - R"), d, n_rep = 5)
svid_perturb <- function(pool, data, gain = 0.03, offset = 0.005, noise = 0.005,
                         soil = 0.02, n_rep = 20L, metric = "auto", seed = 1L) {
  ex <- as.character(pool)
  base <- svid_separability(ex, data, metric = metric)
  set.seed(seed)
  reps <- vapply(seq_len(n_rep), function(r) {
    d2 <- data
    g <- stats::rnorm(1, 1, gain); o <- stats::rnorm(1, 0, offset)
    sb <- stats::rnorm(data$n, 0, soil)
    d2$values <- lapply(data$values, function(m) {
      m * g + o + sb + matrix(stats::rnorm(length(m), 0, noise), nrow(m), ncol(m))
    })
    s <- svid_separability(ex, d2, metric = metric)
    s$score[match(ex, s$expression)]
  }, numeric(length(ex)))
  if (is.null(dim(reps))) reps <- matrix(reps, nrow = 1)
  b <- base$score[match(ex, base$expression)]
  out <- data.frame(expression = ex, score = b,
                    perturbed_mean = rowMeans(reps, na.rm = TRUE),
                    perturbed_sd = apply(reps, 1, stats::sd, na.rm = TRUE),
                    stringsAsFactors = FALSE)
  out$relative_loss <- 1 - out$perturbed_mean / out$score
  out$robustness <- out$perturbed_mean - out$perturbed_sd
  out <- out[order(-out$robustness), ]
  rownames(out) <- NULL
  out
}

#' Transferability across sites, years or scenes
#'
#' Fits on one subset and predicts another, for every pair of levels of a
#' grouping variable. An index that ranks first within a scene but drops
#' sharply when transferred is a local optimum, not a general one, and
#' transferability is usually the property a published index is claimed to have.
#'
#' @param pool A [svid_pool] or character vector.
#' @param data A [svid_data] object.
#' @param by Column of the metadata defining the domains: `"site"`, `"year"`
#'   or `"group"`.
#' @param learner A [svid_learner()].
#' @param n_seeds Seeds per fit.
#'
#' @return A data frame with within-domain and cross-domain scores and the
#'   transfer gap.
#' @export
#' @examples
#' d <- svid_demo_data(sites = 2)
#' svid_transfer(c("(N - R) / (N + R)", "N / (B * R)"), d, by = "site",
#'               n_seeds = 1)
svid_transfer <- function(pool, data, by = c("site", "year", "group"),
                          learner = NULL, n_seeds = 3L) {
  by <- match.arg(by)
  lev <- unique(data$meta[[by]])
  if (length(lev) < 2 || all(is.na(lev))) {
    stop("`", by, "` has fewer than two levels; nothing to transfer between.",
         call. = FALSE)
  }
  if (is.null(learner)) learner <- svid_learner("ranger")
  ex <- as.character(pool)
  vals <- svid_eval(ex, data)
  out <- do.call(rbind, lapply(seq_along(ex), function(i) {
    v <- vals[[i]]
    if (is.null(v)) return(NULL)
    do.call(rbind, lapply(lev, function(l) {
      tr <- which(data$meta[[by]] != l); te <- which(data$meta[[by]] == l)
      if (!length(tr) || !length(te)) return(NULL)
      p <- .fit_predict(learner, v[tr, , drop = FALSE], data$response[tr],
                        v[te, , drop = FALSE], 1L, data$type)
      pin <- .fit_predict(learner, v[tr, , drop = FALSE], data$response[tr],
                          v[tr, , drop = FALSE], 1L, data$type)
      key <- if (data$type == "classification") "oa" else "r2"
      data.frame(expression = ex[i], held_out = l,
                 within = .metrics(data$response[tr], pin, data$type)[[key]],
                 transferred = .metrics(data$response[te], p, data$type)[[key]],
                 stringsAsFactors = FALSE)
    }))
  }))
  out$gap <- out$within - out$transferred
  agg <- stats::aggregate(cbind(within, transferred, gap) ~ expression,
                          data = out, FUN = mean)
  agg <- agg[order(-agg$transferred), ]
  rownames(agg) <- NULL
  structure(agg, class = c("svid_transfer", "data.frame"), detail = out)
}

#' Cross-sensor portability of an index
#'
#' Recomputes an index after resampling the data to the band set of a different
#' platform, so that an index discovered on one sensor can be checked before it
#' is claimed to be general. Bands present on the source sensor but absent on
#' the target are reported, and expressions requiring them are flagged as not
#' portable.
#'
#' @param pool A [svid_pool] or character vector.
#' @param from,to Sensor names or [svid_sensor()] objects.
#'
#' @return A data frame with a `portable` flag and the missing bands.
#' @export
#' @examples
#' \donttest{
#' svid_cross_sensor(c("(N - R) / (N + R)", "(N - S1) / (N + S1)",
#'                     "(RE2 - R) / (RE2 + R)"),
#'                   from = "sentinel2", to = "landsat8")
#' }
svid_cross_sensor <- function(pool, from, to) {
  if (!inherits(from, "svid_sensor")) from <- svid_sensor(from)
  if (!inherits(to, "svid_sensor")) to <- svid_sensor(to)
  ex <- as.character(pool)
  syms <- lapply(ex, function(s) {
    v <- all.vars(str2lang(s))
    v[v %in% .svid_band_table()$symbol]
  })
  miss <- vapply(syms, function(v) paste(setdiff(v, to$symbols), collapse = ","),
                 character(1))
  shift <- vapply(syms, function(v) {
    v <- intersect(v, to$symbols)
    if (!length(v)) return(NA_real_)
    max(abs(from$wavelength[v] - to$wavelength[v]))
  }, numeric(1))
  data.frame(expression = ex, portable = miss == "",
             missing_bands = miss, max_wavelength_shift_nm = shift,
             from = from$name, to = to$name, stringsAsFactors = FALSE)
}

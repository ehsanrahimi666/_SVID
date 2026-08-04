# ---------------------------------------------------------------------------
# Module 9: raster prediction and mapping
# ---------------------------------------------------------------------------

#' Compute an index raster
#'
#' @param cube A [svid_cube()].
#' @param expression A single index expression.
#' @param epoch Epoch index or date. `NULL` returns one layer per epoch.
#' @param filename Optional output path.
#' @return A `SpatRaster`.
#' @export
#' @examples
#' # see vignette("svid-introduction")
#' svid_grammar(c("R", "N"))
svid_index_raster <- function(cube, expression, epoch = NULL, filename = "") {
  .need("terra")
  eps <- if (is.null(epoch)) seq_along(cube$dates) else
    if (is.character(epoch)) match(epoch, cube$dates) else epoch
  nb <- length(cube$bands)
  lyr <- lapply(eps, function(j) {
    env <- new.env(parent = baseenv())
    for (b in seq_len(nb)) {
      assign(cube$bands[b], cube$raster[[(j - 1) * nb + b]] * cube$scale, envir = env)
    }
    eval(str2lang(expression), envir = env)
  })
  out <- terra::rast(lyr)
  names(out) <- cube$dates[eps]
  if (nzchar(filename)) out <- terra::writeRaster(out, filename, overwrite = TRUE)
  out
}

#' Predict a map from a fitted index
#'
#' Trains the learner on the reference samples using one index, then applies it
#' to every pixel of the cube.
#'
#' @param cube A [svid_cube()].
#' @param data A [svid_data] object holding the training samples.
#' @param expression The index expression to use, or `"raw_bands"`.
#' @param learner A [svid_learner()].
#' @param n_seeds Number of seeds; predictions are majority-voted.
#' @param filename Optional output path.
#' @param chunk Rows processed per block.
#'
#' @return A `SpatRaster` of predicted classes or values.
#' @export
#' @examples
#' # see vignette("svid-introduction")
#' svid_learner("ranger")
svid_map <- function(cube, data, expression, learner = NULL, n_seeds = 3L,
                     filename = "", chunk = 200L) {
  .need("terra")
  if (is.null(learner)) learner <- svid_learner("ranger")
  xtr <- if (identical(expression, "raw_bands")) {
    do.call(cbind, lapply(data$bands, function(b) data$values[[b]]))
  } else svid_eval(expression, data)[[1]]
  y <- data$response
  r <- cube$raster
  nb <- length(cube$bands); nd <- length(cube$dates)
  out <- terra::rast(r, nlyrs = 1)
  names(out) <- "prediction"
  bs <- seq(1, terra::nrow(r), by = chunk)
  vals <- vector("list", length(bs))
  for (i in seq_along(bs)) {
    nr <- min(chunk, terra::nrow(r) - bs[i] + 1)
    v <- terra::values(r, row = bs[i], nrows = nr, mat = TRUE) * cube$scale
    env <- new.env(parent = baseenv())
    for (b in seq_len(nb)) {
      assign(cube$bands[b], v[, seq(b, by = nb, length.out = nd), drop = FALSE],
             envir = env)
    }
    xp <- if (identical(expression, "raw_bands")) {
      do.call(cbind, lapply(cube$bands, function(b) get(b, envir = env)))
    } else .clean_feature(matrix(as.numeric(eval(str2lang(expression), envir = env)),
                                nrow(v), nd), "rank")
    ok <- stats::complete.cases(xp)
    p <- rep(NA_real_, nrow(xp))
    if (any(ok)) {
      pr <- lapply(seq_len(n_seeds), function(s)
        .fit_predict(learner, xtr, y, xp[ok, , drop = FALSE], s, data$type))
      p[ok] <- if (data$type == "classification") {
        M <- do.call(cbind, lapply(pr, as.character))
        match(apply(M, 1, function(z) names(which.max(table(z)))), levels(y))
      } else rowMeans(do.call(cbind, pr))
    }
    vals[[i]] <- p
  }
  terra::values(out) <- unlist(vals)
  if (data$type == "classification") {
    levels(out) <- data.frame(id = seq_along(levels(y)), class = levels(y))
  }
  if (nzchar(filename)) out <- terra::writeRaster(out, filename, overwrite = TRUE)
  out
}

#' Area per class in a classified map
#'
#' @param x A classified `SpatRaster`.
#' @param unit `"ha"`, `"km2"` or `"pixels"`.
#' @return A data frame of class areas.
#' @export
#' @examples
#' svid_sensors()
svid_area <- function(x, unit = c("ha", "km2", "pixels")) {
  .need("terra")
  unit <- match.arg(unit)
  f <- terra::freq(x)
  res <- prod(terra::res(x))
  f$area <- switch(unit, pixels = f$count, ha = f$count * res / 1e4,
                   km2 = f$count * res / 1e6)
  f$unit <- unit
  f[, c("value", "count", "area", "unit")]
}

#' Disagreement between two classified maps
#'
#' @param a,b Classified `SpatRaster` objects.
#' @return A list with the disagreement raster, the proportion of pixels that
#'   differ, and a cross-tabulation.
#' @export
#' @examples
#' svid_sensors()
svid_disagree <- function(a, b) {
  .need("terra")
  d <- a != b
  va <- terra::values(a); vb <- terra::values(b)
  ok <- !is.na(va) & !is.na(vb)
  list(raster = d, proportion = mean(va[ok] != vb[ok]),
       table = table(a = va[ok], b = vb[ok]))
}

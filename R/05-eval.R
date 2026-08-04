# ---------------------------------------------------------------------------
# Module 2c: evaluating expressions and building feature matrices
# ---------------------------------------------------------------------------

#' Evaluate index expressions on data
#'
#' @param expressions Character vector of expressions, or a [svid_pool].
#' @param data A [svid_data] object.
#' @param clean How to handle non-finite values. `"rank"` maps `+Inf`, `-Inf`
#'   and `NA` to finite values outside the observed range, preserving rank
#'   order; `"impute"` replaces them with the column mean, reproducing the
#'   common but rank-destroying convention; `"none"` leaves them alone.
#'
#' @return A list of numeric matrices, one per expression, each with samples in
#'   rows and epochs in columns.
#' @export
#' @examples
#' d <- svid_demo_data()
#' v <- svid_eval(c("(N - R) / (N + R)", "N / (B * R)"), d)
#' dim(v[[1]])
#' round(colMeans(v[[1]]), 3)
svid_eval <- function(expressions, data, clean = c("rank", "impute", "none")) {
  clean <- match.arg(clean)
  env <- .band_env(data)
  ex <- as.character(expressions)
  out <- lapply(ex, function(s) {
    v <- tryCatch(eval(str2lang(s), envir = env), error = function(e) NULL)
    if (is.null(v)) return(NULL)
    if (length(v) == 1L) v <- matrix(v, data$n, data$n_epochs)
    v <- matrix(as.numeric(v), data$n, data$n_epochs)
    .clean_feature(v, clean)
  })
  stats::setNames(out, ex)
}

.clean_feature <- function(v, clean = "rank") {
  if (clean == "none") return(v)
  if (clean == "impute") {
    v[!is.finite(v)] <- NA_real_
    mu <- colMeans(v, na.rm = TRUE)
    mu[!is.finite(mu)] <- 0
    idx <- which(is.na(v), arr.ind = TRUE)
    if (nrow(idx)) v[idx] <- mu[idx[, 2]]
    return(pmin(pmax(v, -1e10), 1e10))
  }
  f <- is.finite(v)
  if (!any(f)) return(matrix(0, nrow(v), ncol(v)))
  lo <- min(v[f]); hi <- max(v[f]); sp <- if (hi > lo) hi - lo else 1
  v[is.infinite(v) & v > 0] <- hi + sp
  v[is.infinite(v) & v < 0] <- lo - sp
  v[is.na(v)] <- lo - 2 * sp
  v
}

#' Summarise an index time series into phenological features
#'
#' Reduces a full temporal profile to a small set of interpretable statistics.
#' Useful when acquisition dates are irregular between sites or years, when the
#' number of epochs is large, or when the aim is to compare phenology rather
#' than raw values.
#'
#' @param v A numeric matrix of index values (samples x epochs), or the list
#'   returned by [svid_eval()].
#' @param stats Character vector of statistics to compute: `"mean"`, `"sd"`,
#'   `"min"`, `"max"`, `"amplitude"`, `"integral"`, `"peak_epoch"`,
#'   `"slope_up"`, `"slope_down"`, `"q10"`, `"q90"`.
#'
#' @return A numeric matrix with one column per statistic, or a list of such
#'   matrices.
#' @export
#' @examples
#' d <- svid_demo_data()
#' v <- svid_eval("(N - R) / (N + R)", d)[[1]]
#' head(round(svid_temporal(v), 3), 3)
svid_temporal <- function(v, stats = c("mean", "sd", "min", "max", "amplitude",
                                       "integral", "peak_epoch", "slope_up",
                                       "slope_down")) {
  if (is.list(v)) return(lapply(v, svid_temporal, stats = stats))
  stats <- match.arg(stats, c("mean", "sd", "min", "max", "amplitude",
                              "integral", "peak_epoch", "slope_up",
                              "slope_down", "q10", "q90"), several.ok = TRUE)
  k <- ncol(v)
  pk <- max.col(v, ties.method = "first")
  f <- list(
    mean = rowMeans(v),
    sd = apply(v, 1, stats::sd),
    min = apply(v, 1, min),
    max = apply(v, 1, max),
    amplitude = apply(v, 1, max) - apply(v, 1, min),
    integral = rowSums((v[, -1, drop = FALSE] + v[, -k, drop = FALSE]) / 2),
    peak_epoch = pk,
    slope_up = (v[cbind(seq_len(nrow(v)), pk)] - v[, 1]) / pmax(pk - 1, 1),
    slope_down = (v[, k] - v[cbind(seq_len(nrow(v)), pk)]) / pmax(k - pk, 1),
    q10 = apply(v, 1, stats::quantile, 0.1, names = FALSE),
    q90 = apply(v, 1, stats::quantile, 0.9, names = FALSE)
  )
  out <- do.call(cbind, f[stats])
  colnames(out) <- stats
  out
}

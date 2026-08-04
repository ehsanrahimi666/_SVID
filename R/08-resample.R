# ---------------------------------------------------------------------------
# Module 5a: resampling
# ---------------------------------------------------------------------------

#' Define a resampling scheme
#'
#' Reference pixels digitised inside field parcels are strongly spatially
#' dependent: neighbouring pixels of the same parcel are near-duplicates. A
#' random split then places near-duplicates on both sides of the partition and
#' inflates accuracy. `svid_resample()` therefore offers grouped and spatially
#' blocked schemes, warns when a random split is requested for data that are
#' evidently clustered, and reports the effective number of independent units.
#'
#' @param data A [svid_data] object.
#' @param method One of `"holdout"` (an explicit train/test flag),
#'   `"random"`, `"group"` (leave-groups-out, e.g. parcels), `"spatial_block"`
#'   (contiguous blocks from k-means on coordinates), `"leave_site_out"`,
#'   `"leave_year_out"`.
#' @param folds Number of folds or blocks.
#' @param p Training proportion for `"random"`.
#' @param train Logical vector marking training samples, for `"holdout"`.
#' @param seed Random seed.
#' @param strata Stratify random splits by the response?
#'
#' @return An object of class `svid_resample`, a list of train/test index pairs.
#' @export
#' @examples
#' d <- svid_demo_data()
#' svid_resample(d, "group", folds = 5)
#' svid_resample(d, "spatial_block", folds = 4)
svid_resample <- function(data, method = c("group", "spatial_block", "random",
                                           "holdout", "leave_site_out",
                                           "leave_year_out"),
                          folds = 5L, p = 0.7, train = NULL, seed = 1L,
                          strata = TRUE) {
  method <- match.arg(method)
  n <- data$n
  set.seed(seed)
  grp <- data$meta$group
  ng <- length(unique(grp))
  if (method == "random" && ng < n) {
    warning("the samples fall into only ", ng, " groups but a random split was ",
            "requested; neighbouring samples of the same group will appear in ",
            "both training and test sets. Consider method = \"group\".",
            call. = FALSE)
  }
  splits <- switch(method,
    holdout = {
      if (is.null(train)) stop("`train` is required for method = \"holdout\".",
                               call. = FALSE)
      list(list(train = which(train), test = which(!train)))
    },
    random = {
      idx <- if (strata && data$type == "classification") {
        unlist(lapply(split(seq_len(n), data$response),
                      function(i) sample(i, max(1, floor(p * length(i))))))
      } else sample.int(n, floor(p * n))
      list(list(train = sort(idx), test = setdiff(seq_len(n), idx)))
    },
    group = .fold_by(match(grp, unique(grp)), folds),
    spatial_block = {
      xy <- as.matrix(data$meta[, c("x", "y")])
      if (anyNA(xy)) stop("spatial blocking needs coordinates; supply `coords` ",
                          "when building the data object.", call. = FALSE)
      km <- stats::kmeans(scale(xy), centers = min(folds, nrow(unique(xy))),
                          nstart = 10)
      .fold_by(km$cluster, folds)
    },
    leave_site_out = {
      if (all(is.na(data$meta$site))) stop("no `site` column in the data object.",
                                           call. = FALSE)
      .fold_by(match(data$meta$site, unique(data$meta$site)), Inf)
    },
    leave_year_out = {
      if (all(is.na(data$meta$year))) stop("no `year` column in the data object.",
                                           call. = FALSE)
      .fold_by(match(data$meta$year, unique(data$meta$year)), Inf)
    })
  structure(splits, class = "svid_resample", method = method,
            n = n, n_groups = ng, folds = length(splits))
}

.fold_by <- function(g, folds) {
  u <- unique(g)
  k <- if (is.infinite(folds)) length(u) else min(folds, length(u))
  assign <- stats::setNames(sample(rep_len(seq_len(k), length(u))), u)
  f <- assign[as.character(g)]
  lapply(seq_len(k), function(i) list(train = which(f != i), test = which(f == i)))
}

#' @export
print.svid_resample <- function(x, ...) {
  cat("<svid_resample> method:", attr(x, "method"), "|", length(x), "fold(s)\n")
  cat("  samples:", attr(x, "n"), " independent units:", attr(x, "n_groups"), "\n")
  sz <- vapply(x, function(s) length(s$test), integer(1))
  cat("  test-set sizes:", paste(sz, collapse = ", "), "\n")
  if (attr(x, "n_groups") < attr(x, "n")) {
    cat("  note: effective sample size is the number of independent units,",
        "not the number of pixels\n")
  }
  invisible(x)
}

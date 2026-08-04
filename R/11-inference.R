# ---------------------------------------------------------------------------
# Module 6: inference
# ---------------------------------------------------------------------------

#' Confidence intervals for accuracy estimates
#'
#' Wilson score intervals for proportions (overall and per-class accuracy) and
#' bootstrap intervals for other metrics. Reporting an interval alongside every
#' accuracy makes it immediately visible when a ranking is finer than the data
#' can resolve.
#'
#' @param x A [svid_fit] object.
#' @param level Confidence level.
#' @param metric Metric to bound. Defaults to `"oa"` or `"r2"`.
#' @param n_boot Bootstrap replicates for non-proportion metrics.
#'
#' @return A data frame with `index`, `estimate`, `lower`, `upper`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' fit <- svid_evaluate(c("(N - R) / (N + R)", "N / (B * R)"), d,
#'                      n_seeds = 2, baselines = "none")
#' svid_ci(fit)
svid_ci <- function(x, level = 0.95, metric = NULL, n_boot = 999L) {
  stopifnot(inherits(x, "svid_fit"))
  if (is.null(metric)) metric <- if (x$type == "classification") "oa" else "r2"
  z <- stats::qnorm(1 - (1 - level) / 2)
  m <- x$metrics
  if (metric %in% c("oa", "macro_f1", "balanced_acc") ||
      grepl("^pa_", metric)) {
    n <- m$n_test
    p <- m[[metric]] / 100
    den <- 1 + z^2 / n
    ctr <- p + z^2 / (2 * n)
    mrg <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
    out <- data.frame(index = m$index, estimate = 100 * p,
                      lower = 100 * (ctr - mrg) / den,
                      upper = 100 * (ctr + mrg) / den,
                      metric = metric, stringsAsFactors = FALSE)
  } else {
    tr <- x$truth
    out <- do.call(rbind, lapply(seq_len(nrow(m)), function(i) {
      p <- x$predictions[[m$index[i]]]
      if (is.null(p)) return(NULL)
      b <- replicate(n_boot, {
        s <- sample.int(length(p), replace = TRUE)
        .metrics(tr[s], p[s], x$type)[[metric]]
      })
      data.frame(index = m$index[i], estimate = m[[metric]][i],
                 lower = stats::quantile(b, (1 - level) / 2, names = FALSE),
                 upper = stats::quantile(b, 1 - (1 - level) / 2, names = FALSE),
                 metric = metric, stringsAsFactors = FALSE)
    }))
  }
  rownames(out) <- NULL
  out
}

#' Paired comparison against a reference index
#'
#' McNemar's test for classification and a paired bootstrap for regression,
#' both computed on the same out-of-sample predictions, so that the comparison
#' is paired rather than a difference of two independent estimates.
#'
#' @param x A [svid_fit] object.
#' @param reference Name of the reference index. Defaults to NDVI if present,
#'   otherwise the raw-band baseline.
#' @param method `"exact"` or `"chisq"` for McNemar.
#' @param p_adjust Multiple-comparison adjustment, see [stats::p.adjust()].
#'
#' @return A data frame with the paired test for every index.
#' @export
#' @examples
#' d <- svid_demo_data()
#' fit <- svid_evaluate(c("(N - R) / (N + R)", "N / (B * R)", "N / G"), d,
#'                      n_seeds = 2, baselines = "none")
#' svid_compare(fit, reference = "(N - R) / (N + R)")
svid_compare <- function(x, reference = NULL, method = c("exact", "chisq"),
                         p_adjust = "holm") {
  method <- match.arg(method)
  nm <- names(x$predictions)
  if (is.null(reference)) {
    reference <- if ("(N - R) / (N + R)" %in% nm) "(N - R) / (N + R)" else
      if ("[baseline] NDVI" %in% nm) "[baseline] NDVI" else nm[1]
  }
  if (!reference %in% nm) stop("reference '", reference, "' not among the fitted models.",
                               call. = FALSE)
  tr <- x$truth; ref <- x$predictions[[reference]]
  out <- do.call(rbind, lapply(setdiff(nm, reference), function(k) {
    p <- x$predictions[[k]]
    if (x$type == "classification") {
      a <- p == tr; b <- ref == tr
      n01 <- sum(a & !b); n10 <- sum(!a & b)
      pv <- if (n01 + n10 == 0L) {
        1
      } else if (method == "exact") {
        stats::binom.test(n01, n01 + n10, 0.5)$p.value
      } else {
        s <- (abs(n01 - n10) - 1)^2 / max(n01 + n10, 1)
        stats::pchisq(s, 1, lower.tail = FALSE)
      }
      data.frame(index = k, reference = reference,
                 delta = 100 * (mean(a) - mean(b)),
                 n_better = n01, n_worse = n10, p_value = pv,
                 stringsAsFactors = FALSE)
    } else {
      d <- (ref - tr)^2 - (p - tr)^2
      tt <- stats::t.test(d)
      data.frame(index = k, reference = reference,
                 delta = mean(d), n_better = sum(d > 0), n_worse = sum(d < 0),
                 p_value = tt$p.value, stringsAsFactors = FALSE)
    }
  }))
  out$p_adjusted <- stats::p.adjust(out$p_value, method = p_adjust)
  out <- out[order(-out$delta), ]
  rownames(out) <- NULL
  out
}

#' Indices statistically tied with the best
#'
#' Returns the group of candidates that cannot be separated from the top
#' performer at the requested confidence level, rather than a single winner.
#'
#' @param x A [svid_fit] object.
#' @param alpha Significance level.
#' @param metric Metric used to rank.
#'
#' @return A data frame with a logical `tied` column, of class `svid_tied`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' fit <- svid_evaluate(c("(N - R) / (N + R)", "N / (B * R)", "N / G", "B / G"),
#'                      d, n_seeds = 2, baselines = "none")
#' svid_tied(fit)
svid_tied <- function(x, alpha = 0.05, metric = NULL) {
  ci <- svid_ci(x, level = 1 - alpha, metric = metric)
  best <- which.max(ci$estimate)
  ci$tied <- ci$upper >= ci$lower[best]
  ci <- ci[order(-ci$estimate), ]
  rownames(ci) <- NULL
  structure(ci, class = c("svid_tied", "data.frame"), alpha = alpha)
}

#' @export
print.svid_tied <- function(x, ...) {
  cat("<svid_tied> ", sum(x$tied), " of ", nrow(x),
      " candidates are statistically indistinguishable from the best ",
      "(alpha = ", attr(x, "alpha"), ")\n", sep = "")
  print(utils::head(as.data.frame(x), 10), row.names = FALSE, digits = 4)
  invisible(x)
}

#' Stability of the ranking
#'
#' Quantifies the noise floor of an index ranking from three sources: the
#' spread across random seeds, the spread across resampling folds, and the
#' spread between formulations that belong to the same monotone-equivalence
#' class. The last is the sharpest diagnostic available, because members of one
#' class are mathematically identical and any difference between them is
#' therefore an artefact of the evaluation procedure rather than a property of
#' the index.
#'
#' @param x A [svid_fit] object.
#' @param data A [svid_data] object, required for the within-class component.
#' @param metric Metric to examine.
#' @param max_classes Maximum number of equivalence classes re-evaluated for
#'   the within-class component. Classes are taken in decreasing size, so the
#'   largest and most informative are always included.
#'
#' @return An object of class `svid_stability`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' fit <- svid_evaluate(pool[1:8], d, n_seeds = 3, baselines = "none")
#' svid_stability(fit, d)
svid_stability <- function(x, data = NULL, metric = NULL,
                           max_classes = 25L) {
  if (is.null(metric)) metric <- if (x$type == "classification") "oa" else "r2"
  seed_sd <- x$metrics$sd_seed
  within <- NULL
  if (!is.null(data) && !is.null(x$members)) {
    mem <- x$members[lengths(x$members) > 1]
    if (length(mem) > max_classes) {
      mem <- mem[order(-lengths(mem))][seq_len(max_classes)]
    }
    if (length(mem)) {
      within <- do.call(rbind, lapply(names(mem), function(k) {
        f <- svid_evaluate(mem[[k]], data, resampling = x$resampling,
                           learner = x$learner, n_seeds = 1L, baselines = "none")
        data.frame(class = k, n_members = length(mem[[k]]),
                   spread = diff(range(f$metrics[[metric]])),
                   stringsAsFactors = FALSE)
      }))
    }
  }
  structure(list(metric = metric, seed_sd = seed_sd, within = within,
                 range_top10 = diff(range(utils::head(x$metrics[[metric]], 10)))),
            class = "svid_stability")
}

#' @export
print.svid_stability <- function(x, ...) {
  cat("<svid_stability> metric:", x$metric, "\n")
  cat("  seed-to-seed SD:      median",
      sprintf("%.2f", stats::median(x$seed_sd, na.rm = TRUE)),
      " max", sprintf("%.2f\n", max(x$seed_sd, na.rm = TRUE)))
  cat("  spread across top 10: ", sprintf("%.2f\n", x$range_top10))
  if (!is.null(x$within)) {
    cat("  within-equivalence-class spread (should be 0): median",
        sprintf("%.2f", stats::median(x$within$spread)),
        " max", sprintf("%.2f", max(x$within$spread)),
        "over", nrow(x$within), "classes\n")
    if (max(x$within$spread) >= x$range_top10) {
      cat("  WARNING: identical formulations differ by as much as the whole",
          "top-10 range.\n           The ranking is not resolvable at this",
          "sample size.\n")
    }
  }
  invisible(x)
}

#' Selection optimism
#'
#' Estimates how much the best score in a pool is inflated simply by taking a
#' maximum over many correlated candidates. The estimate is obtained by
#' permuting the response, repeating the whole selection, and recording the
#' best score achievable when no index carries real signal.
#'
#' @param pool A [svid_pool] or character vector.
#' @param data A [svid_data] object.
#' @param n_perm Number of permutations.
#' @param ... Passed to [svid_separability()].
#'
#' @return A list with the observed best score, the null distribution of best
#'   scores, and the implied optimism.
#' @export
#' @examples
#' \donttest{
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' svid_optimism(pool, d, n_perm = 5)
#' }
svid_optimism <- function(pool, data, n_perm = 20L, ...) {
  obs <- max(svid_separability(pool, data, ...)$score, na.rm = TRUE)
  null <- vapply(seq_len(n_perm), function(i) {
    d2 <- data
    set.seed(1000 + i)
    d2$response <- sample(data$response)
    max(svid_separability(pool, d2, ...)$score, na.rm = TRUE)
  }, numeric(1))
  structure(list(observed = obs, null = null,
                 optimism = stats::median(null),
                 p_value = (1 + sum(null >= obs)) / (n_perm + 1),
                 n_candidates = length(pool)),
            class = "svid_optimism")
}

#' @export
print.svid_optimism <- function(x, ...) {
  cat("<svid_optimism> best score over", x$n_candidates, "candidates\n")
  cat("  observed:", sprintf("%.4f", x$observed), "\n")
  cat("  best achievable under a permuted response: median",
      sprintf("%.4f", x$optimism),
      sprintf("[%.4f, %.4f]\n", min(x$null), max(x$null)))
  cat("  permutation p-value:", format.pval(x$p_value), "\n")
  invisible(x)
}

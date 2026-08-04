# ---------------------------------------------------------------------------
# Module 4: fast screening
# ---------------------------------------------------------------------------

#' Class separability of candidate indices
#'
#' Computes cheap statistics that measure how well each index separates the
#' response, without fitting a model. Screening a large pool this way and
#' fitting models only to the best candidates is typically two to three orders
#' of magnitude faster than fitting a model per index, and it is the
#' recommended route for pools of more than a few hundred features.
#'
#' For classification: Jeffries-Matusita distance, transformed divergence,
#' the Fisher ratio and one-versus-rest area under the ROC curve. For
#' regression: absolute Pearson and Spearman correlation with the response.
#'
#' @param pool A [svid_pool] or character vector of expressions.
#' @param data A [svid_data] object.
#' @param metric Statistic to compute. `"auto"` picks `"jm"` for
#'   classification and `"spearman"` for regression.
#' @param aggregate How per-epoch statistics are combined into one score:
#'   `"max"`, `"mean"`, or `"best_epoch"` (which also reports which epoch).
#' @param clean Passed to [svid_eval()].
#'
#' @return A data frame with one row per expression, sorted by decreasing
#'   score, of class `svid_screen_result`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' sep <- svid_separability(pool, d)
#' head(sep, 5)
svid_separability <- function(pool, data,
                              metric = c("auto", "jm", "divergence", "fisher",
                                         "auc", "pearson", "spearman"),
                              aggregate = c("max", "mean", "best_epoch"),
                              clean = "rank") {
  metric <- match.arg(metric); aggregate <- match.arg(aggregate)
  if (metric == "auto") {
    metric <- if (data$type == "classification") "jm" else "spearman"
  }
  ex <- as.character(pool)
  y <- data$response
  vals <- svid_eval(ex, data, clean = clean)
  sc <- vapply(vals, function(v) {
    if (is.null(v)) return(rep(NA_real_, 2))
    s <- vapply(seq_len(ncol(v)), function(j) .sep_one(v[, j], y, metric),
                numeric(1))
    s[!is.finite(s)] <- NA_real_
    if (all(is.na(s))) return(c(NA_real_, NA_real_))
    c(switch(aggregate, max = max(s, na.rm = TRUE),
             mean = mean(s, na.rm = TRUE),
             best_epoch = max(s, na.rm = TRUE)),
      which.max(replace(s, is.na(s), -Inf)))
  }, numeric(2))
  out <- data.frame(expression = ex, score = sc[1, ],
                    best_epoch = data$epochs[sc[2, ]],
                    metric = metric, stringsAsFactors = FALSE)
  out <- out[order(-out$score), ]
  rownames(out) <- NULL
  structure(out, class = c("svid_screen_result", "data.frame"))
}

.sep_one <- function(v, y, metric) {
  if (metric %in% c("pearson", "spearman")) {
    return(abs(suppressWarnings(stats::cor(v, as.numeric(y), method = metric,
                                           use = "complete.obs"))))
  }
  lev <- levels(y); K <- length(lev)
  if (K < 2) return(NA_real_)
  mu <- tapply(v, y, mean); s2 <- tapply(v, y, stats::var)
  s2[!is.finite(s2) | s2 <= 0] <- 1e-12
  pr <- utils::combn(K, 2)
  vals <- vapply(seq_len(ncol(pr)), function(p) {
    i <- pr[1, p]; j <- pr[2, p]
    switch(metric,
      fisher = (mu[i] - mu[j])^2 / (s2[i] + s2[j]),
      {
        b <- 0.125 * (mu[i] - mu[j])^2 * 2 / (s2[i] + s2[j]) +
          0.5 * log((0.5 * (s2[i] + s2[j])) / sqrt(s2[i] * s2[j]))
        if (metric == "jm") 2 * (1 - exp(-b)) else 2 * (1 - exp(-b / 8)) * 1000
      })
  }, numeric(1))
  if (metric == "auc") {
    vals <- vapply(lev, function(l) {
      g <- y == l; r <- rank(v)
      n1 <- sum(g); n0 <- length(v) - n1
      if (n1 == 0 || n0 == 0) return(NA_real_)
      max(c((sum(r[g]) - n1 * (n1 + 1) / 2) / (n1 * n0),
            1 - (sum(r[g]) - n1 * (n1 + 1) / 2) / (n1 * n0)))
    }, numeric(1))
  }
  mean(vals, na.rm = TRUE)
}

#' Take the top candidates from a screening result
#'
#' @param x A [svid_separability()] result.
#' @param n Number of expressions to keep.
#' @param min_score Optional minimum score.
#' @return A [svid_pool].
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' svid_shortlist(svid_separability(pool, d), n = 5)
svid_shortlist <- function(x, n = 50L, min_score = NULL) {
  x <- x[!is.na(x$score), , drop = FALSE]
  if (!is.null(min_score)) x <- x[x$score >= min_score, , drop = FALSE]
  .svid_pool(utils::head(x$expression, n), NULL,
             c(screened = nrow(x), shortlist = min(n, nrow(x))))
}

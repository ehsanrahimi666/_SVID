# ---------------------------------------------------------------------------
# Module 5c: evaluation
# ---------------------------------------------------------------------------

#' Evaluate candidate indices with a model
#'
#' Fits one model per index -- each index contributing its own temporal profile
#' as the feature vector -- under the requested resampling scheme, repeated
#' over several random seeds, and returns pooled out-of-sample predictions and
#' metrics.
#'
#' Two features of the default behaviour are deliberate. Predictions are
#' majority-voted (classification) or averaged (regression) across seeds,
#' because a ranking built on a single random seed is not reproducible: the
#' seed-to-seed spread in overall accuracy is commonly of the same order as the
#' differences separating the leading candidates. And reference baselines are
#' fitted automatically, because an index that beats NDVI but not a model using
#' the raw bands has not demonstrated anything.
#'
#' @param pool A [svid_pool] or character vector of expressions.
#' @param data A [svid_data] object.
#' @param resampling A [svid_resample()] object. Defaults to grouped 5-fold.
#' @param learner A [svid_learner()]. Defaults to `ranger`.
#' @param n_seeds Number of random seeds per index.
#' @param features `"profile"` uses the raw per-epoch values;
#'   `"temporal"` reduces each profile with [svid_temporal()]; `"both"`
#'   concatenates them.
#' @param baselines Character vector of reference models to fit alongside:
#'   any of `"raw_bands"`, `"catalog"`, `"ndvi"`, or `"none"`.
#' @param clean Passed to [svid_eval()].
#' @param progress Print progress?
#'
#' @return An object of class `svid_fit`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' fit <- svid_evaluate(pool[1:6], d, n_seeds = 2, baselines = "raw_bands")
#' fit
#' head(as.data.frame(fit), 4)
svid_evaluate <- function(pool, data, resampling = NULL, learner = NULL,
                          n_seeds = 5L, features = c("profile", "temporal", "both"),
                          baselines = c("raw_bands"), clean = "rank",
                          progress = FALSE) {
  features <- match.arg(features)
  if (is.null(learner)) learner <- svid_learner("ranger")
  if (is.null(resampling)) resampling <- svid_resample(data, "group", folds = 5L)
  ex <- as.character(pool)
  seeds <- seq_len(n_seeds) * 17L + 25L

  build <- function(v) {
    switch(features, profile = v, temporal = svid_temporal(v),
           both = cbind(v, svid_temporal(v)))
  }
  vals <- svid_eval(ex, data, clean = clean)
  feats <- lapply(vals, function(v) if (is.null(v)) NULL else build(v))

  extra <- list()
  if (!identical(baselines, "none")) {
    if ("raw_bands" %in% baselines) {
      extra[["[baseline] all bands"]] <-
        do.call(cbind, lapply(data$bands, function(b) build(data$values[[b]])))
    }
    if ("ndvi" %in% baselines && all(c("N", "R") %in% data$bands)) {
      extra[["[baseline] NDVI"]] <- build(svid_eval("(N - R) / (N + R)", data, clean)[[1]])
    }
    if ("catalog" %in% baselines) {
      cg <- svid_catalog(data, as_pool = TRUE)
      cv <- svid_eval(cg, data, clean = clean)
      nm <- attr(cg, "source")
      for (i in seq_along(cv)) {
        if (!is.null(cv[[i]]))
          extra[[paste0("[catalog] ", nm[i])]] <- build(cv[[i]])
      }
    }
  }
  feats <- c(feats, extra)

  res <- vector("list", length(feats))
  preds <- vector("list", length(feats))
  for (i in seq_along(feats)) {
    if (is.null(feats[[i]])) next
    pm <- .pooled_predict(feats[[i]], data, resampling, learner, seeds)
    preds[[i]] <- pm$vote
    res[[i]] <- cbind(data.frame(index = names(feats)[i], stringsAsFactors = FALSE),
                      .metrics(pm$truth, pm$vote, data$type),
                      data.frame(sd_seed = pm$sd_seed, n_test = length(pm$truth)))
    if (progress && i %% 25 == 0) message("  evaluated ", i, "/", length(feats))
  }
  out <- do.call(rbind, res)
  key <- if (data$type == "classification") "kappa" else "r2"
  out <- out[order(-out[[key]]), ]
  out$rank <- seq_len(nrow(out))
  rownames(out) <- NULL
  names(preds) <- names(feats)
  structure(list(metrics = out, predictions = preds[!vapply(preds, is.null, logical(1))],
                 truth = data$response, type = data$type, learner = learner,
                 resampling = resampling, n_seeds = n_seeds,
                 members = attr(pool, "members"), features = features),
            class = "svid_fit")
}

.pooled_predict <- function(x, data, resampling, learner, seeds) {
  y <- data$response
  per_seed <- lapply(seeds, function(s) {
    p <- if (is.factor(y)) factor(rep(NA, data$n), levels = levels(y)) else rep(NA_real_, data$n)
    for (sp in resampling) {
      p[sp$test] <- .fit_predict(learner, x[sp$train, , drop = FALSE], y[sp$train],
                                 x[sp$test, , drop = FALSE], s, data$type)
    }
    p
  })
  used <- sort(unique(unlist(lapply(resampling, `[[`, "test"))))
  if (is.factor(y)) {
    M <- do.call(cbind, lapply(per_seed, as.character))
    vote <- factor(apply(M[used, , drop = FALSE], 1, function(z) {
      t <- table(z); names(t)[which.max(t)]
    }), levels = levels(y))
    acc <- vapply(per_seed, function(p) mean(p[used] == y[used], na.rm = TRUE), numeric(1))
  } else {
    M <- do.call(cbind, per_seed)
    vote <- rowMeans(M[used, , drop = FALSE], na.rm = TRUE)
    acc <- vapply(per_seed, function(p) stats::cor(p[used], y[used], use = "complete.obs"), numeric(1))
  }
  list(vote = vote, truth = y[used], idx = used,
       sd_seed = if (length(acc) > 1) stats::sd(acc) * 100 else NA_real_)
}

.metrics <- function(truth, pred, type) {
  if (type == "classification") {
    lev <- levels(truth)
    cm <- table(factor(truth, lev), factor(pred, lev))
    oa <- sum(diag(cm)) / sum(cm)
    pe <- sum(rowSums(cm) * colSums(cm)) / sum(cm)^2
    pa <- diag(cm) / pmax(rowSums(cm), 1)
    ua <- diag(cm) / pmax(colSums(cm), 1)
    f1 <- 2 * pa * ua / pmax(pa + ua, 1e-12)
    out <- data.frame(oa = 100 * oa, kappa = (oa - pe) / (1 - pe),
                      macro_f1 = 100 * mean(f1, na.rm = TRUE),
                      balanced_acc = 100 * mean(pa, na.rm = TRUE))
    pa_df <- as.data.frame(t(100 * pa))
    names(pa_df) <- paste0("pa_", lev)
    cbind(out, pa_df)
  } else {
    e <- pred - truth
    data.frame(rmse = sqrt(mean(e^2)), mae = mean(abs(e)),
               r2 = 1 - sum(e^2) / sum((truth - mean(truth))^2),
               bias = mean(e))
  }
}

#' Fit reference baselines only
#'
#' @inheritParams svid_evaluate
#' @param which Which baselines to fit.
#' @param ... Further arguments passed to [svid_evaluate()].
#' @return A [svid_fit] object containing only the baselines.
#' @export
#' @examples
#' d <- svid_demo_data()
#' svid_baseline(d, which = c("raw_bands", "ndvi"), n_seeds = 2)
svid_baseline <- function(data, which = c("raw_bands", "ndvi", "catalog"),
                          resampling = NULL, learner = NULL, n_seeds = 5L, ...) {
  svid_evaluate(character(0), data, resampling = resampling, learner = learner,
                n_seeds = n_seeds, baselines = which, ...)
}

#' Evaluate candidate indices against a continuous response
#'
#' Convenience wrapper around [svid_evaluate()] for regression targets such as
#' leaf area index, canopy chlorophyll, biomass or yield. Reports RMSE, MAE,
#' bias and R-squared.
#'
#' @inheritParams svid_evaluate
#' @param ... Further arguments passed to [svid_evaluate()].
#' @return An object of class `svid_fit`.
#' @export
#' @examples
#' sim <- svid_simulate(n = 120, sensor = "sentinel2", response = "lai", seed = 1)
#' pool <- svid_prune(svid_generate(svid_grammar(c("R", "N", "S1"))), sim)
#' svid_regress(pool[1:5], sim, n_seeds = 2, baselines = "none")
svid_regress <- function(pool, data, ...) {
  if (data$type != "regression") {
    stop("`data` has a categorical response; use svid_evaluate() instead.",
         call. = FALSE)
  }
  svid_evaluate(pool, data, ...)
}

#' @export
print.svid_fit <- function(x, ...) {
  cat("<svid_fit>", x$type, "|", nrow(x$metrics), "models |",
      x$n_seeds, "seeds |", attr(x$resampling, "method"), "resampling\n")
  key <- if (x$type == "classification") c("oa", "kappa", "macro_f1") else c("r2", "rmse")
  m <- x$metrics[, c("rank", "index", key, "sd_seed"), drop = FALSE]
  print(utils::head(m, 8), row.names = FALSE, digits = 4)
  if (nrow(m) > 8) cat("  ... ", nrow(m) - 8, " more\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.svid_fit <- function(x, ...) x$metrics

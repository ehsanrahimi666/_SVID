# ---------------------------------------------------------------------------
# Module 5b: learners
# ---------------------------------------------------------------------------

#' Define the learner used for evaluation
#'
#' A thin, pluggable adapter. `svid` does not implement learning algorithms; it
#' wraps existing ones so that the same evaluation protocol can be applied to
#' any of them. The choice matters for redundancy reduction: rank-based
#' learners such as random forests are invariant to monotone transformations,
#' which is what [svid_equivalence()] exploits, whereas linear models are not.
#'
#' @param type One of `"ranger"`, `"randomForest"`, `"rpart"`, `"lm"`,
#'   `"glm"`, or `"custom"`.
#' @param ... Arguments passed to the underlying fitting function, e.g.
#'   `num.trees`, `mtry`, `max.depth` for `ranger`.
#' @param fit,predict For `type = "custom"`, functions with signatures
#'   `fit(x, y, ...)` and `predict(model, newx)`.
#' @param invariance The transformation invariance of the learner, used as the
#'   default for [svid_prune()]. Inferred from `type` if `NULL`.
#'
#' @return An object of class `svid_learner`.
#' @export
#' @examples
#' svid_learner("ranger", num.trees = 300, max.depth = 10)
#' svid_learner("lm")
svid_learner <- function(type = c("ranger", "randomForest", "rpart", "lm",
                                  "glm", "custom"),
                         ..., fit = NULL, predict = NULL, invariance = NULL) {
  type <- match.arg(type)
  if (is.null(invariance)) {
    invariance <- if (type %in% c("ranger", "randomForest", "rpart")) "rank" else "affine"
  }
  structure(list(type = type, args = list(...), fit = fit, predict = predict,
                 invariance = invariance), class = "svid_learner")
}

#' @export
print.svid_learner <- function(x, ...) {
  cat("<svid_learner>", x$type, "| invariance:", x$invariance, "\n")
  if (length(x$args)) {
    cat("  args:", paste(names(x$args), unlist(lapply(x$args, format)),
                         sep = "=", collapse = " "), "\n")
  }
  invisible(x)
}

.fit_predict <- function(learner, xtr, ytr, xte, seed = 1L, type) {
  set.seed(seed)
  colnames(xtr) <- colnames(xte) <- paste0("f", seq_len(ncol(xtr)))
  switch(learner$type,
    ranger = {
      .need("ranger")
      args <- utils::modifyList(list(num.trees = 500, min.node.size = 1,
                                     seed = seed, num.threads = 1),
                                learner$args)
      m <- do.call(ranger::ranger,
                   c(list(x = as.data.frame(xtr), y = ytr,
                          classification = identical(type, "classification")),
                     args))
      stats::predict(m, data = as.data.frame(xte), num.threads = 1)$predictions
    },
    randomForest = {
      .need("randomForest")
      m <- do.call(randomForest::randomForest,
                   c(list(x = as.data.frame(xtr), y = ytr), learner$args))
      stats::predict(m, as.data.frame(xte))
    },
    rpart = {
      .need("rpart")
      df <- data.frame(y = ytr, xtr)
      m <- do.call(rpart::rpart, c(list(formula = y ~ ., data = df),
                                   learner$args))
      p <- stats::predict(m, as.data.frame(xte),
                          type = if (type == "classification") "class" else "vector")
      p
    },
    lm = {
      df <- data.frame(y = as.numeric(ytr), xtr)
      m <- stats::lm(y ~ ., data = df)
      stats::predict(m, as.data.frame(xte))
    },
    glm = {
      df <- data.frame(y = ytr, xtr)
      m <- do.call(stats::glm, c(list(formula = y ~ ., data = df), learner$args))
      p <- stats::predict(m, as.data.frame(xte), type = "response")
      if (type == "classification") factor(levels(ytr)[1 + (p > 0.5)], levels(levels(ytr))) else p
    },
    custom = {
      m <- learner$fit(xtr, ytr)
      learner$predict(m, xte)
    })
}

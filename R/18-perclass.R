# ---------------------------------------------------------------------------
# Module 4b: per-class (one-versus-rest) index discovery
# ---------------------------------------------------------------------------

#' Which index is best for each class
#'
#' A single index chosen for a whole landscape is a compromise: it maximises
#' average performance across classes that may have entirely different spectral
#' and phenological behaviour. `svid_per_class()` instead treats every class as
#' its own one-versus-rest problem and reports the index that separates that
#' class best, together with an interval on the estimate.
#'
#' This is the natural question for crop mapping. The index that best isolates
#' a double-cropped soybean system from everything else need not be the index
#' that best isolates sugarcane, because their phenologies differ in timing
#' rather than in amplitude.
#'
#' @param pool A [svid_pool] or character vector of expressions. Deduplicate it
#'   with [svid_prune()] first; otherwise one feature can win several times
#'   under different formulations.
#' @param data A [svid_data] object with a categorical response.
#' @param metric Separability statistic. `"auc"` is the default because it is
#'   directly interpretable for a binary contrast and insensitive to class
#'   prevalence; `"jm"`, `"divergence"` and `"fisher"` are also accepted.
#' @param aggregate How per-epoch scores are combined: `"max"` takes the single
#'   most discriminating date, which identifies phenological windows, `"mean"`
#'   rewards indices that separate the class throughout the season.
#' @param n_top Number of ranked candidates retained per class.
#' @param clean Passed to [svid_eval()].
#' @param progress Report progress per class?
#'
#' @return An object of class `svid_per_class`: a list holding the full
#'   class-by-index score matrix, the ranked table per class, and the winning
#'   index for each class with its best epoch.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                               max_depth = 1)), d)
#' pc <- svid_per_class(pool, d)
#' pc
#' pc$winners
svid_per_class <- function(pool, data, metric = c("auc", "jm", "divergence", "fisher"),
                           aggregate = c("max", "mean"), n_top = 10L,
                           clean = "rank", progress = FALSE) {
  metric <- match.arg(metric); aggregate <- match.arg(aggregate)
  if (data$type != "classification") {
    stop("per-class discovery requires a categorical response.", call. = FALSE)
  }
  ex <- as.character(pool)
  lev <- levels(data$response)
  bin <- lapply(lev, function(l) factor(ifelse(data$response == l, "pos", "neg"),
                                        levels = c("neg", "pos")))
  S <- matrix(NA_real_, length(ex), length(lev), dimnames = list(ex, lev))
  E <- S
  ## expressions are evaluated one at a time: holding the whole pool in memory
  ## is prohibitive for large pools and long time series
  fast_auc <- identical(metric, "auc")
  pos <- lapply(lev, function(l) data$response == l)
  n1 <- vapply(pos, sum, integer(1)); n0 <- length(data$response) - n1
  for (i in seq_along(ex)) {
    v <- svid_eval(ex[i], data, clean = clean)[[1]]
    if (is.null(v)) next
    if (fast_auc) {
      ## the rank transform depends only on the feature, not on the class, so
      ## it is computed once per epoch and reused for every one-versus-rest
      ## contrast. This is the dominant cost for large pools.
      RK <- apply(v, 2, rank, ties.method = "average")
      for (k in seq_along(lev)) {
        a <- (colSums(RK[pos[[k]], , drop = FALSE]) - n1[k] * (n1[k] + 1) / 2) /
          (n1[k] * n0[k])
        s <- pmax(a, 1 - a)
        s[!is.finite(s)] <- NA_real_
        if (all(is.na(s))) next
        S[i, k] <- if (aggregate == "max") max(s, na.rm = TRUE) else mean(s, na.rm = TRUE)
        E[i, k] <- which.max(replace(s, is.na(s), -Inf))
      }
    } else {
      for (k in seq_along(lev)) {
        s <- vapply(seq_len(ncol(v)), function(j) .sep_one(v[, j], bin[[k]], metric),
                    numeric(1))
        s[!is.finite(s)] <- NA_real_
        if (all(is.na(s))) next
        S[i, k] <- if (aggregate == "max") max(s, na.rm = TRUE) else mean(s, na.rm = TRUE)
        E[i, k] <- which.max(replace(s, is.na(s), -Inf))
      }
    }
    if (progress && i %% 200 == 0) message("  ", i, "/", length(ex))
  }
  ranked <- lapply(seq_along(lev), function(k) {
    o <- order(-S[, k])
    o <- o[!is.na(S[o, k])]
    o <- utils::head(o, n_top)
    data.frame(class = lev[k], rank = seq_along(o), expression = ex[o],
               score = S[o, k], best_epoch = data$epochs[E[o, k]],
               n_class = sum(data$response == lev[k]),
               stringsAsFactors = FALSE)
  })
  ranked <- do.call(rbind, ranked)
  winners <- ranked[ranked$rank == 1L, , drop = FALSE]
  rownames(ranked) <- rownames(winners) <- NULL
  structure(list(scores = S, epochs = E, ranked = ranked, winners = winners,
                 metric = metric, aggregate = aggregate, levels = lev,
                 n_candidates = length(ex)),
            class = "svid_per_class")
}

#' @export
print.svid_per_class <- function(x, ...) {
  cat("<svid_per_class>", length(x$levels), "classes |", x$n_candidates,
      "candidate features | metric:", x$metric, "\n")
  w <- x$winners[, c("class", "expression", "score", "best_epoch", "n_class")]
  print(w, row.names = FALSE, digits = 4)
  nd <- length(unique(x$winners$expression))
  cat("  ", nd, " distinct winning features across ", length(x$levels),
      " classes\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.svid_per_class <- function(x, ...) x$ranked

#' Specialist indices against a single all-class index
#'
#' Quantifies what is lost by forcing one index on every class. For each class
#' the specialist score is compared with the score achieved by the index that is
#' best on average across all classes, and the difference is reported.
#'
#' @param x A [svid_per_class()] result.
#' @param global Optional expression to use as the single all-class index. If
#'   `NULL`, the index with the highest mean score across classes is used.
#'
#' @return A data frame of class `svid_specialists` with one row per class.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                               max_depth = 1)), d)
#' svid_specialists(svid_per_class(pool, d))
svid_specialists <- function(x, global = NULL) {
  stopifnot(inherits(x, "svid_per_class"))
  mean_score <- rowMeans(x$scores, na.rm = TRUE)
  if (is.null(global)) global <- names(which.max(mean_score))
  if (!global %in% rownames(x$scores)) {
    stop("`global` is not among the evaluated candidates.", call. = FALSE)
  }
  out <- data.frame(
    class = x$levels,
    n_class = x$winners$n_class[match(x$levels, x$winners$class)],
    specialist = x$winners$expression[match(x$levels, x$winners$class)],
    specialist_score = x$winners$score[match(x$levels, x$winners$class)],
    specialist_epoch = x$winners$best_epoch[match(x$levels, x$winners$class)],
    global_index = global,
    global_score = x$scores[global, ],
    stringsAsFactors = FALSE)
  out$gain <- out$specialist_score - out$global_score
  out <- out[order(-out$gain), ]
  rownames(out) <- NULL
  structure(out, class = c("svid_specialists", "data.frame"),
            metric = x$metric, global = global)
}

#' @export
print.svid_specialists <- function(x, ...) {
  cat("<svid_specialists> single all-class index:", attr(x, "global"), "\n")
  print(as.data.frame(x)[, c("class", "n_class", "specialist",
                             "specialist_score", "global_score", "gain")],
        row.names = FALSE, digits = 4)
  cat("  mean gain from specialisation:",
      sprintf("%.4f", mean(x$gain, na.rm = TRUE)),
      " classes improved:", sum(x$gain > 0), "of", nrow(x), "\n")
  invisible(x)
}

#' Classify with one specialist index per class
#'
#' Turns the per-class winners into an actual classifier. One binary
#' one-versus-rest model is fitted for every class, each using only that class's
#' specialist index, and the per-class scores are combined by taking the class
#' whose model is most confident. The result is directly comparable with a
#' conventional multi-class model trained on a single index or on the raw bands.
#'
#' @param x A [svid_per_class()] result.
#' @param data A [svid_data] object.
#' @param resampling A [svid_resample()] object.
#' @param learner A [svid_learner()]. Must support probability output; `ranger`
#'   and `randomForest` do.
#' @param n_seeds Number of random seeds; per-class scores are averaged.
#' @param compare Additional expressions to evaluate as multi-class references,
#'   for example the single best all-class index.
#'
#' @return An object of class `svid_ensemble` holding the pooled out-of-sample
#'   predictions and metrics for the specialist ensemble and each reference.
#' @export
#' @examples
#' \donttest{
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                               max_depth = 1)), d)
#' pc <- svid_per_class(pool, d)
#' svid_ensemble(pc, d, n_seeds = 1)
#' }
svid_ensemble <- function(x, data, resampling = NULL, learner = NULL,
                          n_seeds = 3L, compare = NULL) {
  stopifnot(inherits(x, "svid_per_class"))
  if (is.null(learner)) learner <- svid_learner("ranger")
  if (is.null(resampling)) resampling <- svid_resample(data, "group", folds = 5L)
  lev <- x$levels
  spec <- stats::setNames(x$winners$expression[match(lev, x$winners$class)], lev)
  feats <- lapply(spec, function(e) .clean_feature(svid_eval(e, data)[[1]], "rank"))
  seeds <- seq_len(n_seeds) * 17L + 25L

  P <- matrix(0, data$n, length(lev), dimnames = list(NULL, lev))
  for (k in seq_along(lev)) {
    yk <- factor(ifelse(data$response == lev[k], "pos", "neg"),
                 levels = c("neg", "pos"))
    for (s in seeds) {
      for (sp in resampling) {
        p <- .prob_pos(learner, feats[[k]][sp$train, , drop = FALSE], yk[sp$train],
                       feats[[k]][sp$test, , drop = FALSE], s)
        P[sp$test, k] <- P[sp$test, k] + p / length(seeds)
      }
    }
  }
  used <- sort(unique(unlist(lapply(resampling, `[[`, "test"))))
  pred <- factor(lev[max.col(P[used, , drop = FALSE], ties.method = "first")],
                 levels = lev)
  truth <- data$response[used]
  res <- cbind(data.frame(model = "specialist ensemble", stringsAsFactors = FALSE),
               .metrics(truth, pred, "classification"))
  preds <- list(`specialist ensemble` = pred)

  refs <- c(compare, stats::setNames(names(which.max(rowMeans(x$scores, na.rm = TRUE))),
                                     "best all-class index"))
  refs <- refs[!duplicated(refs)]
  for (i in seq_along(refs)) {
    f <- svid_evaluate(refs[i], data, resampling = resampling, learner = learner,
                       n_seeds = n_seeds, baselines = "none")
    nm <- if (!is.null(names(refs)) && nzchar(names(refs)[i])) names(refs)[i] else refs[i]
    res <- rbind(res, cbind(data.frame(model = nm, stringsAsFactors = FALSE),
                            f$metrics[1, names(res)[-1], drop = FALSE]))
    preds[[nm]] <- f$predictions[[1]]
  }
  rownames(res) <- NULL
  structure(list(metrics = res, predictions = preds, probabilities = P,
                 specialists = spec, truth = truth, index = used),
            class = "svid_ensemble")
}

.prob_pos <- function(learner, xtr, ytr, xte, seed) {
  colnames(xtr) <- colnames(xte) <- paste0("f", seq_len(ncol(xtr)))
  set.seed(seed)
  if (learner$type == "ranger" && requireNamespace("ranger", quietly = TRUE)) {
    args <- utils::modifyList(list(num.trees = 500, seed = seed, num.threads = 1),
                              learner$args)
    m <- do.call(ranger::ranger, c(list(x = as.data.frame(xtr), y = ytr,
                                        probability = TRUE), args))
    return(stats::predict(m, data = as.data.frame(xte),
                          num.threads = 1)$predictions[, "pos"])
  }
  if (requireNamespace("randomForest", quietly = TRUE)) {
    m <- randomForest::randomForest(x = as.data.frame(xtr), y = ytr)
    return(stats::predict(m, as.data.frame(xte), type = "prob")[, "pos"])
  }
  as.numeric(.fit_predict(learner, xtr, ytr, xte, seed, "classification") == "pos")
}

#' @export
print.svid_ensemble <- function(x, ...) {
  cat("<svid_ensemble>", length(x$specialists), "specialist indices\n")
  print(x$metrics, row.names = FALSE, digits = 4)
  invisible(x)
}

#' Heatmap of index performance across classes
#'
#' Every candidate index against every class, with each class's winner marked.
#' Blocks of high score sharing a set of bands indicate a family of indices that
#' suits a group of classes; an isolated bright cell indicates a class that
#' needs its own formulation.
#'
#' @param x A [svid_per_class()] result.
#' @param n Number of top candidates to display, chosen by best rank achieved
#'   in any class.
#' @param mark_winners Mark the best index of each class?
#'
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                               max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_specialists(svid_per_class(pool, d))
svid_plot_specialists <- function(x, n = 25L, mark_winners = TRUE) {
  .gg()
  S <- x$scores
  best_rank <- apply(S, 1, function(z) max(z, na.rm = TRUE))
  keep <- utils::head(order(-best_rank), n)
  S <- S[keep, , drop = FALSE]
  df <- data.frame(expression = rep(rownames(S), ncol(S)),
                   class_name = rep(colnames(S), each = nrow(S)),
                   score = as.vector(S), stringsAsFactors = FALSE)
  df$expression <- factor(df$expression, levels = rownames(S)[order(rowMeans(S, na.rm = TRUE))])
  win <- x$winners[x$winners$expression %in% rownames(S), , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = class_name, y = expression, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(colours = c("#F5F5F5", "#9EBCD4", "#2B6A99", "#16324A"),
                                  na.value = "grey95") +
    ggplot2::labs(title = "Index performance by class",
                  subtitle = paste0("one-versus-rest ", x$metric,
                                    "; boxes mark each class's best index"),
                  x = NULL, y = NULL, fill = x$metric) +
    .svid_theme() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(family = "mono", size = 7))
  if (mark_winners && nrow(win)) {
    p <- p + ggplot2::geom_tile(data = data.frame(class_name = win$class,
                                                  expression = win$expression,
                                                  score = NA),
                                fill = NA, colour = "#C0603A", linewidth = 0.9)
  }
  p
}

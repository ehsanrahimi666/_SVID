# ---------------------------------------------------------------------------
# Module 8: guided search
# ---------------------------------------------------------------------------

#' Guided search over index space
#'
#' Exhaustive enumeration becomes impractical beyond about six bands or depth
#' three. `svid_search()` explores the same space with an evolutionary
#' algorithm guided by a separability objective and a parsimony penalty,
#' canonicalising every candidate so that the population cannot fill up with
#' rewritings of one feature.
#'
#' @param data A [svid_data] object.
#' @param bands Band symbols to combine. Defaults to those in `data`.
#' @param n_generations,population Search effort.
#' @param max_depth Maximum expression depth.
#' @param metric Objective, passed to [svid_separability()].
#' @param parsimony Penalty per band occurrence, subtracted from the score.
#' @param elitism Proportion of the population carried forward unchanged.
#' @param mutation Probability of mutating a subtree.
#' @param seed Random seed.
#' @param verbose Report the best score each generation?
#'
#' @return An object of class `svid_search` containing the final population,
#'   the best expressions and the trace of best scores.
#' @export
#' @examples
#' d <- svid_demo_data()
#' s <- svid_search(d, n_generations = 5, population = 40, seed = 1)
#' s
#' head(s$best, 3)
svid_search <- function(data, bands = NULL, n_generations = 20L,
                        population = 200L, max_depth = 3L, metric = "auto",
                        parsimony = 0.002, elitism = 0.1, mutation = 0.3,
                        seed = 1L, verbose = FALSE) {
  if (is.null(bands)) bands <- data$bands
  set.seed(seed)
  ops <- c("+", "-", "*", "/")
  rand_expr <- function(depth) {
    if (depth <= 1L || stats::runif(1) < 0.3) return(sample(bands, 1))
    sprintf("(%s) %s (%s)", rand_expr(depth - 1L), sample(ops, 1),
            rand_expr(depth - 1L))
  }
  pop <- unique(vapply(seq_len(population), function(i)
    rand_expr(sample.int(max_depth, 1)), character(1)))
  n_elite <- max(2L, floor(elitism * population))
  trace <- numeric(n_generations)
  for (g in seq_len(n_generations)) {
    pop <- unique(pop[!is.na(svid_canonicalize(pop))])
    sc <- svid_separability(pop, data, metric = metric)
    nt <- vapply(sc$expression, function(s) .terms_in(str2lang(s)), integer(1))
    sc$fitness <- sc$score - parsimony * nt
    sc <- sc[order(-sc$fitness), ]
    trace[g] <- sc$fitness[1]
    if (verbose) message("gen ", g, ": best ", signif(sc$fitness[1], 4),
                         "  ", sc$expression[1])
    elite <- utils::head(sc$expression, n_elite)
    kids <- character(0)
    while (length(kids) < population - n_elite) {
      a <- sample(utils::head(sc$expression, max(4, nrow(sc) %/% 2)), 1)
      if (stats::runif(1) < mutation) {
        kids <- c(kids, sprintf("(%s) %s (%s)", a, sample(ops, 1),
                                rand_expr(max(1L, max_depth - 1L))))
      } else {
        b <- sample(sc$expression, 1)
        kids <- c(kids, sprintf("(%s) %s (%s)", a, sample(ops, 1), b))
      }
    }
    pop <- unique(c(elite, kids))
    pop <- pop[vapply(pop, function(s)
      .terms_in(str2lang(s)) <= 2L * max_depth, logical(1))]
  }
  final <- svid_separability(pop, data, metric = metric)
  structure(list(best = final, population = pop, trace = trace,
                 n_generations = n_generations, bands = bands),
            class = "svid_search")
}

#' @export
print.svid_search <- function(x, ...) {
  cat("<svid_search>", x$n_generations, "generations |",
      length(x$population), "final candidates\n")
  cat("  best score:", signif(max(x$best$score, na.rm = TRUE), 4), "\n")
  print(utils::head(x$best[, c("expression", "score")], 5), row.names = FALSE)
  invisible(x)
}

#' Beam search over index space
#'
#' A deterministic alternative to [svid_search()]: keeps the `width` best
#' expressions at each depth and extends only those. Faster and reproducible,
#' at the cost of exploring less.
#'
#' @param data A [svid_data] object.
#' @param bands Band symbols. Defaults to those in `data`.
#' @param width Beam width.
#' @param max_depth Maximum depth.
#' @param metric Objective, passed to [svid_separability()].
#'
#' @return A [svid_separability()] result for the final beam.
#' @export
#' @examples
#' \donttest{
#' d <- svid_demo_data()
#' head(svid_beam(d, width = 10, max_depth = 2), 5)
#' }
svid_beam <- function(data, bands = NULL, width = 30L, max_depth = 3L,
                      metric = "auto") {
  if (is.null(bands)) bands <- data$bands
  beam <- bands
  for (d in seq_len(max_depth)) {
    grid <- expand.grid(a = beam, op = c("+", "-", "*", "/"), b = bands,
                        stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    cand <- unique(c(beam, sprintf("(%s) %s (%s)", grid$a, grid$op, grid$b)))
    cand <- cand[!duplicated(svid_canonicalize(cand))]
    sc <- svid_separability(cand, data, metric = metric)
    beam <- utils::head(sc$expression, width)
  }
  svid_separability(beam, data, metric = metric)
}

#' Pareto front of accuracy against complexity
#'
#' Returns the non-dominated set: expressions that no other expression beats on
#' both accuracy and simplicity. Reporting a front rather than a single winner
#' is the honest presentation when the leading candidates are statistically
#' tied.
#'
#' @param x A [svid_fit] object or a [svid_separability()] result.
#' @param metric Metric to maximise. Inferred if `NULL`.
#'
#' @return A data frame with a logical `front` column.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' svid_pareto(svid_separability(pool, d))
svid_pareto <- function(x, metric = NULL) {
  if (inherits(x, "svid_fit")) {
    if (is.null(metric)) metric <- if (x$type == "classification") "oa" else "r2"
    df <- data.frame(expression = x$metrics$index, score = x$metrics[[metric]],
                     stringsAsFactors = FALSE)
  } else {
    df <- data.frame(expression = x$expression, score = x$score,
                     stringsAsFactors = FALSE)
  }
  df <- df[!is.na(df$score), , drop = FALSE]
  df$complexity <- vapply(df$expression, function(s)
    tryCatch(.terms_in(str2lang(s)), error = function(e) NA_integer_),
    integer(1), USE.NAMES = FALSE)
  df <- df[!is.na(df$complexity), , drop = FALSE]
  df$front <- vapply(seq_len(nrow(df)), function(i) {
    !any(df$score >= df$score[i] & df$complexity <= df$complexity[i] &
           (df$score > df$score[i] | df$complexity < df$complexity[i]))
  }, logical(1))
  df <- df[order(df$complexity, -df$score), ]
  rownames(df) <- NULL
  df
}

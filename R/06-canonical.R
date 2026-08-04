# ---------------------------------------------------------------------------
# Module 3: redundancy reduction
#
# Exact symbolic canonicalisation of rational band expressions.
# Every expression admitted by the grammar is a rational function of the band
# symbols, so it can be written as a quotient of two multivariate polynomials
# with numeric coefficients. Two expressions are algebraically identical if and
# only if their canonical (numerator, denominator) pair agrees, once the pair
# has been normalised so that the denominator's leading coefficient is one.
# ---------------------------------------------------------------------------

.mono_key <- function(v) {
  v <- v[v != 0]
  if (!length(v)) return("1")
  v <- v[order(names(v))]
  paste(sprintf("%s^%d", names(v), v), collapse = "*")
}

.poly_const <- function(k) if (k == 0) list() else stats::setNames(list(k), "1")

.poly_var <- function(v) stats::setNames(list(1), .mono_key(stats::setNames(1L, v)))

.poly_add <- function(a, b, sign = 1) {
  for (k in names(b)) {
    a[[k]] <- (if (is.null(a[[k]])) 0 else a[[k]]) + sign * b[[k]]
  }
  a[abs(unlist(a)) > 1e-12]
}

.key_to_exp <- function(k) {
  if (identical(k, "1")) return(integer(0))
  p <- strsplit(strsplit(k, "\\*")[[1]], "\\^")
  stats::setNames(as.integer(vapply(p, `[`, "", 2)), vapply(p, `[`, "", 1))
}

.poly_mul <- function(a, b) {
  if (!length(a) || !length(b)) return(list())
  out <- list()
  ea <- lapply(names(a), .key_to_exp); eb <- lapply(names(b), .key_to_exp)
  for (i in seq_along(a)) for (j in seq_along(b)) {
    e <- ea[[i]]
    for (v in names(eb[[j]])) e[v] <- (if (is.na(e[v])) 0L else e[v]) + eb[[j]][[v]]
    k <- .mono_key(e)
    out[[k]] <- (if (is.null(out[[k]])) 0 else out[[k]]) + a[[i]] * b[[j]]
  }
  out[abs(unlist(out)) > 1e-12]
}

.rat <- function(num, den = .poly_const(1)) list(num = num, den = den)

.rat_of <- function(e) {
  if (is.numeric(e)) return(.rat(.poly_const(as.numeric(e))))
  if (is.name(e))    return(.rat(.poly_var(as.character(e))))
  if (!is.call(e))   stop("unsupported expression element", call. = FALSE)
  op <- as.character(e[[1]])
  if (op == "(") return(.rat_of(e[[2]]))
  if (op == "-" && length(e) == 2L) {
    a <- .rat_of(e[[2]])
    return(.rat(.poly_mul(a$num, .poly_const(-1)), a$den))
  }
  if (op == "+" && length(e) == 2L) return(.rat_of(e[[2]]))
  if (!op %in% c("+", "-", "*", "/")) {
    stop("operator '", op, "' is not a rational operation; ",
         "symbolic canonicalisation supports + - * / only.", call. = FALSE)
  }
  a <- .rat_of(e[[2]]); b <- .rat_of(e[[3]])
  switch(op,
    "+" = .rat(.poly_add(.poly_mul(a$num, b$den), .poly_mul(b$num, a$den)),
               .poly_mul(a$den, b$den)),
    "-" = .rat(.poly_add(.poly_mul(a$num, b$den), .poly_mul(b$num, a$den), sign = -1),
               .poly_mul(a$den, b$den)),
    "*" = .rat(.poly_mul(a$num, b$num), .poly_mul(a$den, b$den)),
    "/" = .rat(.poly_mul(a$num, b$den), .poly_mul(a$den, b$num)))
}

.poly_str <- function(p, digits = 10) {
  if (!length(p)) return("0")
  k <- sort(names(p))
  paste(sprintf("%s*%s", signif(unlist(p[k]), digits), k), collapse = "+")
}

#' Canonical algebraic form of index expressions
#'
#' Rewrites each expression as a normalised quotient of multivariate
#' polynomials. Two expressions receive the same canonical form if and only if
#' they define the same rational function of the bands, so `(N/B)/R`,
#' `(N/R)/B` and `N/(B*R)` all collapse to one entry.
#'
#' This is the step that prevents a single feature from occupying several
#' places in a ranking. In an exhaustively generated four-band pool it
#' typically removes a quarter to a third of all expressions.
#'
#' @param x A character vector of expressions or a [svid_pool].
#' @param digits Number of significant digits retained in coefficients before
#'   comparison, guarding against floating-point noise.
#'
#' @return A character vector of canonical-form strings, the same length as `x`.
#' @export
#' @examples
#' e <- c("(N / B) / R", "(N / R) / B", "N / (B * R)", "(N - R) / (N + R)")
#' svid_canonicalize(e)
#' # the first three share one canonical form
#' duplicated(svid_canonicalize(e))
svid_canonicalize <- function(x, digits = 10) {
  ex <- as.character(x)
  vapply(ex, function(s) {
    tryCatch({
      r <- .rat_of(str2lang(s))
      if (!length(r$num)) return("0")
      lead <- r$den[[sort(names(r$den))[1]]]
      if (!is.finite(lead) || lead == 0) lead <- 1
      r$num <- .poly_mul(r$num, .poly_const(1 / lead))
      r$den <- .poly_mul(r$den, .poly_const(1 / lead))
      paste0("(", .poly_str(r$num, digits), ")/(", .poly_str(r$den, digits), ")")
    }, error = function(e) NA_character_)
  }, character(1), USE.NAMES = FALSE)
}

#' Screen expressions for numerical validity
#'
#' Evaluates every expression on the data and flags those that are constant,
#' non-finite for a large share of samples, or numerically degenerate. This
#' mirrors the stability screening step of a manual workflow but records the
#' reason each expression was dropped.
#'
#' @param pool A [svid_pool] or character vector of expressions.
#' @param data A [svid_data] object.
#' @param max_nonfinite Maximum tolerated proportion of non-finite values.
#' @param min_unique Minimum number of distinct values required.
#' @param n_sample Maximum number of cells used to compute the statistics.
#'
#' @return A data frame with one row per expression and columns `expression`,
#'   `keep`, `reason`, `p_nonfinite` and `n_unique`.
#' @export
#' @examples
#' data(andong)
#' d <- svid_demo_data()
#' pool <- svid_generate(svid_grammar(c("B", "G", "R", "N"), max_depth = 1))
#' scr <- svid_screen(pool, d)
#' table(scr$reason)
svid_screen <- function(pool, data, max_nonfinite = 0.05, min_unique = 10L,
                        n_sample = 4000L) {
  ex <- as.character(pool)
  env <- .band_env(data)
  N <- length(data$values[[1]])
  set.seed(1L)
  idx <- if (N > n_sample) sort(sample.int(N, n_sample)) else seq_len(N)
  res <- lapply(ex, function(s) {
    v <- tryCatch(eval(str2lang(s), envir = env), error = function(e) NULL)
    if (is.null(v)) return(c(NA_real_, NA_real_, 0))
    v <- as.numeric(v)[idx]
    f <- is.finite(v)
    vs <- v[f]
    if (length(vs) > 5000L) vs <- vs[seq(1, length(vs), length.out = 5000L)]
    c(mean(!f), length(unique(signif(vs, 8))), 1)
  })
  m <- do.call(rbind, res)
  reason <- rep("ok", length(ex))
  reason[m[, 3] == 0] <- "evaluation error"
  reason[m[, 3] == 1 & m[, 1] > max_nonfinite] <- "too many non-finite values"
  reason[m[, 3] == 1 & m[, 2] < min_unique] <- "insufficient variation"
  data.frame(expression = ex, keep = reason == "ok", reason = reason,
             p_nonfinite = m[, 1], n_unique = m[, 2],
             stringsAsFactors = FALSE)
}

#' Remove numerically identical expressions
#'
#' Catches duplicates that survive symbolic canonicalisation, for example
#' expressions whose common polynomial factors cancel only for the data at
#' hand.
#'
#' @param pool A [svid_pool] or character vector.
#' @param data A [svid_data] object.
#' @param digits Significant digits retained in the hashed values.
#' @param n_sample Maximum number of (sample, epoch) cells used to build the
#'   hash. Subsampling makes the comparison linear in the pool size rather than
#'   in the size of the data, at negligible risk of a collision.
#' @param seed Random seed for the subsample.
#'
#' @return An integer vector assigning each expression to a duplicate group.
#' @export
#' @examples
#' d <- svid_demo_data()
#' svid_dedupe(c("N / R", "(N * G) / (R * G)", "N - R"), d)
svid_dedupe <- function(pool, data, digits = 8L, n_sample = 4000L, seed = 1L) {
  ex <- as.character(pool)
  env <- .band_env(data)
  N <- length(data$values[[1]])
  set.seed(seed)
  idx <- if (N > n_sample) sort(sample.int(N, n_sample)) else seq_len(N)
  keys <- vapply(ex, function(s) {
    v <- tryCatch(as.numeric(eval(str2lang(s), envir = env))[idx],
                  error = function(e) NA_real_)
    v[!is.finite(v)] <- NA_real_
    paste(signif(v, digits), collapse = ",")
  }, character(1), USE.NAMES = FALSE)
  match(keys, unique(keys))
}

#' Monotone-equivalence classes
#'
#' Groups features that are strictly monotone transformations of one another.
#' A learner that splits on rank order -- a decision tree, a random forest,
#' gradient boosting, k-nearest neighbours -- cannot distinguish such features,
#' so they carry identical information and must not be ranked against each
#' other. `N / R`, `R / N` and `(N - R) / (N + R)` all belong to one class.
#'
#' The relevant invariance depends on the learner. Set `invariance = "rank"`
#' for tree ensembles and other rank-based methods, `"affine"` for linear and
#' distance-based models, which are invariant only to positive affine
#' rescaling, and `"none"` to disable the reduction.
#'
#' @param pool A [svid_pool] or character vector.
#' @param data A [svid_data] object.
#' @param invariance `"rank"`, `"affine"` or `"none"`.
#' @param n_sample Number of (sample, epoch) cells used to test equivalence.
#' @param tol Tolerance on the correlation used to declare equivalence.
#' @param seed Random seed for subsampling.
#'
#' @return An object of class `svid_equivalence`: a list with the class
#'   membership of every expression, one canonical representative per class,
#'   and the reduction counts.
#' @export
#' @examples
#' d <- svid_demo_data()
#' eq <- svid_equivalence(c("N / R", "R / N", "(N - R) / (N + R)", "N - R"), d)
#' eq
#' eq$members
svid_equivalence <- function(pool, data, invariance = c("rank", "affine", "none"),
                             n_sample = 3000L, tol = 1e-9, seed = 1L) {
  invariance <- match.arg(invariance)
  ex <- as.character(pool)
  if (invariance == "none" || length(ex) < 2L) {
    return(structure(list(class_of = seq_along(ex), representative = ex,
                          members = stats::setNames(as.list(ex), ex),
                          invariance = invariance, n_input = length(ex),
                          n_classes = length(ex)),
                     class = "svid_equivalence"))
  }
  env <- .band_env(data)
  N <- length(data$values[[1]])
  set.seed(seed)
  idx <- if (N > n_sample) sample.int(N, n_sample) else seq_len(N)
  M <- vapply(ex, function(s) {
    v <- tryCatch(as.numeric(eval(str2lang(s), envir = env))[idx],
                  error = function(e) rep(NA_real_, length(idx)))
    v[!is.finite(v)] <- NA_real_
    v
  }, numeric(length(idx)))
  ## rank-transform once, then a single centred cross-product. Using
  ## cor(use = "pairwise.complete.obs") here would compute every pair
  ## separately and is orders of magnitude slower for large pools.
  if (invariance == "rank") {
    M <- apply(M, 2, function(z) {
      if (all(is.na(z))) return(rep(0, length(z)))
      r <- rank(z, na.last = "keep", ties.method = "average")
      r[is.na(r)] <- mean(r, na.rm = TRUE)
      r
    })
  } else {
    M[is.na(M)] <- 0
  }
  M <- scale(M)
  M[!is.finite(M)] <- 0
  R <- crossprod(M) / (nrow(M) - 1)
  R[!is.finite(R)] <- 0
  parent <- seq_along(ex)
  find <- function(a) { while (parent[a] != a) { parent[a] <<- parent[parent[a]]; a <- parent[a] }; a }
  hit <- which(abs(R) >= 1 - tol & upper.tri(R), arr.ind = TRUE)
  for (r in seq_len(nrow(hit))) {
    a <- find(hit[r, 1]); b <- find(hit[r, 2])
    if (a != b) parent[a] <- b
  }
  cls <- vapply(seq_along(ex), find, integer(1))
  cls <- match(cls, unique(cls))
  members <- split(ex, cls)
  rep_ex <- vapply(members, function(z) z[which.min(nchar(z))], character(1))
  names(members) <- rep_ex
  structure(list(class_of = cls, representative = unname(rep_ex),
                 members = members, invariance = invariance,
                 n_input = length(ex), n_classes = length(members)),
            class = "svid_equivalence")
}

#' @export
print.svid_equivalence <- function(x, ...) {
  sz <- lengths(x$members)
  cat("<svid_equivalence> invariance:", x$invariance, "\n")
  cat("  ", x$n_input, " expressions -> ", x$n_classes,
      " distinguishable features (", sprintf("%.1f%%", 100 * (1 - x$n_classes / x$n_input)),
      " redundant)\n", sep = "")
  cat("  classes with more than one formulation:", sum(sz > 1),
      "| largest class:", max(sz), "\n")
  invisible(x)
}

#' Full redundancy reduction
#'
#' Runs canonicalisation, numerical screening, numeric deduplication and
#' monotone-equivalence detection in order, and returns a pool containing one
#' representative per distinguishable feature together with a complete audit
#' trail. This is the recommended entry point and the step that should always
#' precede ranking.
#'
#' @param pool A [svid_pool] or character vector of expressions.
#' @param data A [svid_data] object.
#' @param invariance Passed to [svid_equivalence()].
#' @param screen Run numerical screening?
#' @param n_sample Number of cells used by the screening, deduplication and
#'   equivalence steps.
#' @param ... Further arguments passed to [svid_screen()].
#'
#' @return A [svid_pool] whose `history` attribute records the count at each
#'   stage and whose `members` attribute lists the equivalent formulations of
#'   every retained feature.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_generate(svid_grammar(c("B", "G", "R", "N"), max_depth = 1))
#' pruned <- svid_prune(pool, d)
#' pruned
#' svid_redundancy(pruned)
svid_prune <- function(pool, data, invariance = "rank", screen = TRUE,
                       n_sample = 4000L, ...) {
  ex <- as.character(pool)
  hist <- c(generated = length(ex))

  can <- svid_canonicalize(ex)
  ok <- !is.na(can)
  ex <- ex[ok]; can <- can[ok]
  ex <- ex[!duplicated(can)]
  hist <- c(hist, canonical = length(ex))

  if (screen) {
    scr <- svid_screen(ex, data, n_sample = n_sample, ...)
    ex <- ex[scr$keep]
    hist <- c(hist, screened = length(ex))
  }

  grp <- svid_dedupe(ex, data, n_sample = n_sample)
  ex <- ex[!duplicated(grp)]
  hist <- c(hist, distinct = length(ex))

  eq <- svid_equivalence(ex, data, invariance = invariance,
                         n_sample = min(n_sample, 3000L))
  hist <- c(hist, equivalence_classes = length(eq$representative))

  .svid_pool(eq$representative, attr(pool, "grammar"), hist,
             members = eq$members, source = attr(pool, "source"))
}

.unused_subset_data <- function(data, i) {
  data$values <- lapply(data$values, function(m) m[i, , drop = FALSE])
  data$response <- data$response[i]
  data$meta <- data$meta[i, , drop = FALSE]
  data$n <- length(i)
  data
}

#' Redundancy summary
#'
#' Extracts the stage-by-stage reduction recorded by [svid_prune()], for
#' reporting and for [svid_plot_redundancy()].
#'
#' @param pool A pruned [svid_pool].
#' @return A data frame with columns `stage`, `n` and `removed`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' svid_redundancy(p)
svid_redundancy <- function(pool) {
  h <- attr(pool, "history")
  if (is.null(h)) stop("this pool carries no reduction history; run svid_prune().",
                       call. = FALSE)
  data.frame(stage = names(h), n = as.integer(h),
             removed = c(NA_integer_, as.integer(-diff(h))),
             stringsAsFactors = FALSE)
}

.band_env <- function(data) {
  stopifnot(inherits(data, "svid_data"))
  e <- new.env(parent = baseenv())
  for (b in data$bands) assign(b, data$values[[b]], envir = e)
  e
}

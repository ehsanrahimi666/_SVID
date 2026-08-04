# ---------------------------------------------------------------------------
# Module 2: expression generation
# ---------------------------------------------------------------------------

#' Define the search grammar
#'
#' Controls the shape of the candidate expression pool: which bands may appear,
#' which arithmetic operators are allowed, how deep expressions may nest, and
#' whether numeric coefficients are admitted. The size of the pool grows very
#' quickly with `max_depth` and with the number of bands, so
#' [svid_grammar_size()] is provided to check before generating.
#'
#' @param bands Character vector of band symbols to combine.
#' @param ops Allowed binary operators, any of `"+"`, `"-"`, `"*"`, `"/"`.
#' @param max_depth Maximum nesting depth. Depth 1 gives `a op b`; depth 2
#'   allows one operand to be itself a depth-1 expression; and so on.
#' @param coefficients Optional numeric vector of constants that may multiply a
#'   band, as in `N - 0.66 * S1`. Keep this short: each value multiplies the
#'   pool size.
#' @param allow_self Allow a band to appear more than once in an expression?
#' @param normalised_only Restrict the pool to normalised-difference forms
#'   `(a - b) / (a + b)` over all band pairs. This reproduces the classical
#'   two-band NDSI search of the hyperspectral literature and is much smaller
#'   than the free grammar.
#' @param max_terms Maximum number of band occurrences in an expression.
#'
#' @return An object of class `svid_grammar`.
#' @export
#' @examples
#' g <- svid_grammar(bands = c("B", "G", "R", "N"), max_depth = 2)
#' g
#' svid_grammar_size(g)
#'
#' # the classical two-band normalised-difference search, on a SWIR sensor
#' g2 <- svid_grammar(bands = svid_bands(sensor = "landsat8")$symbol,
#'                    normalised_only = TRUE)
#' svid_grammar_size(g2)
svid_grammar <- function(bands, ops = c("+", "-", "*", "/"), max_depth = 2L,
                         coefficients = NULL, allow_self = FALSE,
                         normalised_only = FALSE, max_terms = 4L) {
  ops <- match.arg(ops, c("+", "-", "*", "/"), several.ok = TRUE)
  bad <- setdiff(bands, .svid_band_table()$symbol)
  if (length(bad)) stop("unknown band symbol(s): ", paste(bad, collapse = ", "),
                        call. = FALSE)
  if (max_depth < 1L) stop("`max_depth` must be at least 1.", call. = FALSE)
  structure(list(bands = bands, ops = ops, max_depth = as.integer(max_depth),
                 coefficients = coefficients, allow_self = allow_self,
                 normalised_only = normalised_only,
                 max_terms = as.integer(max_terms)),
            class = "svid_grammar")
}

#' Size of the pool a grammar will generate
#'
#' Counts the expressions without building them, so an over-ambitious grammar
#' can be caught before it exhausts memory.
#'
#' @param grammar A [svid_grammar()].
#' @return The number of expressions, or `Inf` when the count overflows.
#' @export
#' @examples
#' svid_grammar_size(svid_grammar(c("B", "G", "R", "N"), max_depth = 2))
svid_grammar_size <- function(grammar) {
  if (grammar$normalised_only) {
    nb <- length(grammar$bands)
    return(as.integer(nb * (nb - 1) / 2))
  }
  nb <- length(grammar$bands) * (1 + length(grammar$coefficients))
  no <- length(grammar$ops)
  cnt <- nb
  for (d in seq_len(grammar$max_depth)) {
    cnt <- nb + no * cnt * cnt
    if (!is.finite(cnt) || cnt > 1e15) return(Inf)
  }
  cnt
}

#' @export
print.svid_grammar <- function(x, ...) {
  cat("<svid_grammar>\n")
  cat("  bands:", paste(x$bands, collapse = " "), "\n")
  if (x$normalised_only) {
    cat("  form: normalised difference (a - b)/(a + b) over all band pairs\n")
  } else {
    cat("  operators:", paste(x$ops, collapse = " "),
        " max depth:", x$max_depth, " max terms:", x$max_terms, "\n")
    if (length(x$coefficients))
      cat("  coefficients:", paste(x$coefficients, collapse = ", "), "\n")
  }
  sz <- svid_grammar_size(x)
  cat("  upper bound on pool size:",
      if (is.finite(sz)) format(sz, big.mark = ",") else "overflow", "\n")
  invisible(x)
}

.terms_in <- function(e) {
  if (is.name(e)) return(1L)
  if (is.numeric(e)) return(0L)
  sum(vapply(as.list(e)[-1], .terms_in, integer(1)))
}

#' Generate a pool of candidate index expressions
#'
#' Enumerates every expression admitted by `grammar`, up to the requested
#' depth. The result is a character vector of R expressions in which band
#' symbols appear as variables, so any expression can be evaluated directly by
#' [svid_eval()] or pasted into other software.
#'
#' Generation is deliberately separated from redundancy reduction: the raw pool
#' is expected to contain duplicates, and [svid_prune()] removes them while
#' recording how many were removed at each stage.
#'
#' @param grammar A [svid_grammar()].
#' @param max_size Safety limit. Generation stops with an error if the pool
#'   would exceed this many expressions; raise it deliberately.
#' @param verbose Print progress?
#'
#' @return An object of class `svid_pool`: a character vector of expressions
#'   with attributes recording the grammar and generation history.
#' @export
#' @examples
#' pool <- svid_generate(svid_grammar(c("B", "G", "R", "N"), max_depth = 1))
#' pool
#' head(as.character(pool), 8)
#'
#' # a two-band normalised-difference pool including SWIR
#' nd <- svid_generate(svid_grammar(c("B", "G", "R", "N", "S1", "S2"),
#'                                  normalised_only = TRUE))
#' as.character(nd)
svid_generate <- function(grammar, max_size = 2e5, verbose = FALSE) {
  stopifnot(inherits(grammar, "svid_grammar"))
  if (grammar$normalised_only) {
    cb <- utils::combn(grammar$bands, 2)
    ex <- sprintf("(%s - %s) / (%s + %s)", cb[1, ], cb[2, ], cb[1, ], cb[2, ])
    return(.svid_pool(ex, grammar, c(generated = length(ex))))
  }
  sz <- svid_grammar_size(grammar)
  if (!is.finite(sz) || sz > max_size) {
    stop("this grammar would generate about ",
         if (is.finite(sz)) format(sz, big.mark = ",") else "an astronomical number of",
         " expressions, above `max_size` = ", format(max_size, big.mark = ","),
         ".\nReduce `max_depth`, use fewer bands, set `normalised_only = TRUE`, ",
         "or use svid_search() for guided search instead of enumeration.",
         call. = FALSE)
  }
  atoms <- grammar$bands
  if (length(grammar$coefficients)) {
    atoms <- c(atoms, as.vector(outer(grammar$coefficients, grammar$bands,
                                      function(k, b) paste0(k, " * ", b))))
  }
  level <- atoms
  all_ex <- character(0)
  for (d in seq_len(grammar$max_depth)) {
    grid <- expand.grid(a = level, op = grammar$ops, b = level,
                        stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
    new <- sprintf("(%s) %s (%s)", grid$a, grid$op, grid$b)
    if (!grammar$allow_self) {
      bl <- lapply(strsplit(gsub("[^A-Za-z0-9 ]", " ", new), " +"), function(z) {
        z <- z[z %in% grammar$bands]; anyDuplicated(z) == 0L
      })
      new <- new[unlist(bl)]
    }
    keep <- vapply(new, function(s) {
      tryCatch(.terms_in(str2lang(s)) <= grammar$max_terms,
               error = function(e) FALSE)
    }, logical(1), USE.NAMES = FALSE)
    new <- new[keep]
    all_ex <- unique(c(all_ex, new))
    level <- unique(c(atoms, new))
    if (verbose) message("depth ", d, ": ", length(all_ex), " expressions")
    if (length(all_ex) > max_size) {
      stop("pool exceeded `max_size` at depth ", d, ".", call. = FALSE)
    }
  }
  .svid_pool(all_ex, grammar, c(generated = length(all_ex)))
}

.svid_pool <- function(expressions, grammar = NULL, history = NULL,
                       classes = NULL, members = NULL, source = NULL) {
  structure(unname(expressions), class = "svid_pool", grammar = grammar,
            history = history, classes = classes, members = members,
            source = source)
}

#' @export
print.svid_pool <- function(x, ...) {
  cat("<svid_pool>", length(x), "expressions\n")
  h <- attr(x, "history")
  if (!is.null(h) && length(h) > 1) {
    cat("  reduction:", paste(sprintf("%s=%d", names(h), h), collapse = " -> "), "\n")
  }
  m <- attr(x, "members")
  if (!is.null(m)) {
    sz <- lengths(m)
    cat("  equivalence classes:", length(m),
        sprintf("(largest %d formulations, %d classes with >1)",
                max(sz), sum(sz > 1)), "\n")
  }
  cat("  e.g. ", paste(utils::head(as.character(x), 3), collapse = " | "), "\n", sep = "")
  invisible(x)
}

#' @export
as.character.svid_pool <- function(x, ...) as.vector(unclass(x))

#' @export
`[.svid_pool` <- function(x, i, ...) {
  y <- unclass(x)[i]
  m <- attr(x, "members")
  .svid_pool(y, attr(x, "grammar"), attr(x, "history"),
             members = if (is.null(m)) NULL else m[y], source = attr(x, "source"))
}

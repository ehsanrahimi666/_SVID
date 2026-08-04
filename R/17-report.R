# ---------------------------------------------------------------------------
# Module 10b: reporting and export
# ---------------------------------------------------------------------------

#' Summarise a discovery run
#'
#' Collects the pool reduction, the ranking, the tied set, the stability
#' diagnostics and the baseline comparison into a single object suitable for
#' printing, for `knitr` tables, or for [svid_report()].
#'
#' @param pool A pruned [svid_pool].
#' @param fit A [svid_fit] object.
#' @param data A [svid_data] object.
#' @param alpha Significance level for the tied set.
#'
#' @return An object of class `svid_result`.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' fit <- svid_evaluate(pool[1:5], d, n_seeds = 2, baselines = "raw_bands")
#' svid_result(pool, fit, d)
svid_result <- function(pool, fit, data, alpha = 0.05) {
  structure(list(
    redundancy = svid_redundancy(pool),
    ranking = fit$metrics,
    tied = svid_tied(fit, alpha = alpha),
    stability = svid_stability(fit, NULL),
    comparison = tryCatch(svid_compare(fit), error = function(e) NULL),
    data = utils::capture.output(print(data)),
    resampling = utils::capture.output(print(fit$resampling)),
    n_features = length(pool)
  ), class = "svid_result")
}

#' @export
print.svid_result <- function(x, ...) {
  cat("<svid_result>\n\n-- pool reduction --\n")
  print(x$redundancy, row.names = FALSE)
  cat("\n-- ranking (top 8) --\n")
  print(utils::head(x$ranking, 8), row.names = FALSE, digits = 4)
  cat("\n-- statistically tied with the best --\n")
  cat("  ", sum(x$tied$tied), "of", nrow(x$tied), "candidates\n")
  cat("\n-- stability --\n"); print(x$stability)
  invisible(x)
}

#' Render an HTML report
#'
#' Writes a self-contained report with every diagnostic and figure produced by
#' a discovery run, so that the analysis behind a published index can be
#' distributed alongside it.
#'
#' @param result A [svid_result()] object.
#' @param file Output path.
#' @param title Report title.
#' @param pool,fit,data Optional objects used to regenerate figures.
#' @return The output path, invisibly.
#' @export
#' @examples
#' # writes to a temporary file
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' fit <- svid_evaluate(pool[1:4], d, n_seeds = 1, baselines = "none")
#' res <- svid_result(pool, fit, d)
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   f <- svid_report(res, file.path(tempdir(), "svid.html"),
#'                    pool = pool, fit = fit, data = d)
#' }
svid_report <- function(result, file = "svid-report.html",
                        title = "Spectral index discovery report",
                        pool = NULL, fit = NULL, data = NULL) {
  .need("rmarkdown")
  rmd <- tempfile(fileext = ".Rmd")
  writeLines(c(
    "---", paste0("title: \"", title, "\""),
    "output: html_document", "---", "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,",
    "                      fig.width = 8, fig.height = 4.5, dpi = 130)",
    "library(svid)", "```", "",
    "## Pool reduction", "```{r}", "knitr::kable(result$redundancy)", "```",
    "```{r}", "if (!is.null(pool)) svid_plot_redundancy(pool)", "```", "",
    "## Ranking", "```{r}", "knitr::kable(head(result$ranking, 15), digits = 3)", "```",
    "```{r}", "if (!is.null(fit)) svid_plot_ranking(fit)", "```", "",
    "## Statistical ties", "```{r}",
    "knitr::kable(head(as.data.frame(result$tied), 15), digits = 3)", "```", "",
    "## Stability", "```{r}", "result$stability", "```", "",
    "## Per-class accuracy", "```{r}",
    "if (!is.null(fit) && fit$type == 'classification') svid_plot_perclass(fit)",
    "```", "",
    "## Paired comparisons", "```{r}",
    "if (!is.null(result$comparison)) knitr::kable(result$comparison, digits = 4)",
    "```", "",
    "## Data and resampling", "```{r}", "cat(result$data, sep = '\\n')",
    "cat(result$resampling, sep = '\\n')", "```"), rmd)
  e <- new.env()
  assign("result", result, envir = e); assign("pool", pool, envir = e)
  assign("fit", fit, envir = e); assign("data", data, envir = e)
  rmarkdown::render(rmd, output_file = normalizePath(file, mustWork = FALSE),
                    envir = e, quiet = TRUE)
  invisible(file)
}

#' Export discovered indices
#'
#' Writes the retained expressions in a portable form, optionally in the syntax
#' of another platform so that a discovered index can be applied outside R.
#'
#' @param pool A [svid_pool], or a [svid_fit] whose ranking is exported.
#' @param file Output path. `NULL` returns the text.
#' @param format `"csv"`, `"json"`, `"gee"` (Earth Engine JavaScript) or
#'   `"python"`.
#' @param sensor Optional sensor used to translate band symbols to native names.
#' @param n Number of expressions to export.
#'
#' @return The written text, invisibly.
#' @export
#' @examples
#' \donttest{
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' cat(svid_export(pool[1:3], NULL, format = "gee", sensor = "sentinel2"))
#' }
svid_export <- function(pool, file = NULL, format = c("csv", "json", "gee", "python"),
                        sensor = NULL, n = Inf) {
  format <- match.arg(format)
  ex <- if (inherits(pool, "svid_fit")) pool$metrics$index else as.character(pool)
  ex <- utils::head(ex, n)
  nat <- function(s) {
    if (is.null(sensor)) return(s)
    sn <- if (inherits(sensor, "svid_sensor")) sensor else svid_sensor(sensor)
    for (b in names(sort(nchar(sn$symbols), decreasing = TRUE))) {}
    for (b in sn$symbols[order(-nchar(sn$symbols))]) {
      s <- gsub(paste0("\\b", b, "\\b"), sn$native[[b]], s)
    }
    s
  }
  txt <- switch(format,
    csv = paste(c("expression", paste0("\"", ex, "\"")), collapse = "\n"),
    json = paste0("[\n", paste(sprintf("  {\"expression\": \"%s\"}", ex),
                               collapse = ",\n"), "\n]"),
    gee = paste(c("// Generated by the svid R package",
                  sprintf("var svid_%02d = function(img) { return img.expression('%s', {%s}); };",
                          seq_along(ex), nat(ex),
                          vapply(ex, function(s) {
                            v <- intersect(all.vars(str2lang(s)), .svid_band_table()$symbol)
                            paste(sprintf("'%s': img.select('%s')", nat(v), nat(v)),
                                  collapse = ", ")
                          }, character(1)))), collapse = "\n"),
    python = paste(c("# Generated by the svid R package", "import numpy as np",
                     sprintf("def svid_%02d(%s):\n    return %s", seq_along(ex),
                             vapply(ex, function(s)
                               paste(intersect(all.vars(str2lang(s)),
                                               .svid_band_table()$symbol),
                                     collapse = ", "), character(1)),
                             ex)), collapse = "\n\n"))
  if (!is.null(file)) writeLines(txt, file)
  invisible(txt)
}

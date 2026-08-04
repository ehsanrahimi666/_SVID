# ---------------------------------------------------------------------------
# Module 2b: published index catalogue
# ---------------------------------------------------------------------------

#' Published spectral indices
#'
#' Returns entries from the bundled Awesome Spectral Indices catalogue
#' (Montero et al. 2023), filtered to the bands actually available. Use it to
#' benchmark discovered indices against the published literature rather than
#' against NDVI alone: a discovered index that beats NDVI but not the best of
#' 280 existing indices is not a discovery.
#'
#' @param bands Character vector of available band symbols, or a [svid_data]
#'   or [svid_sensor] object from which they are taken.
#' @param domain Optional application domains to keep, e.g. `"vegetation"`,
#'   `"water"`, `"soil"`, `"burn"`, `"urban"`, `"snow"`.
#' @param names Optional character vector of short names to select directly.
#' @param as_pool Return a [svid_pool] rather than a data frame?
#'
#' @return A data frame of catalogue entries, or a [svid_pool].
#' @source Montero, D., Aybar, C., Mahecha, M.D. et al. (2023) A standardized
#'   catalogue of spectral indices to advance the use of remote sensing in
#'   Earth system research. \emph{Scientific Data} 10, 197.
#'   \doi{10.1038/s41597-023-02096-0}
#' @export
#' @examples
#' head(svid_catalog(c("B", "G", "R", "N")), 4)
#' nrow(svid_catalog(c("B", "G", "R", "N", "S1", "S2"), domain = "vegetation"))
#' svid_catalog(names = c("NDVI", "EVI", "NDWI"))
svid_catalog <- function(bands = NULL, domain = NULL, names = NULL,
                         as_pool = FALSE) {
  cat_df <- get("asi_indices", envir = asNamespace("svid"))
  if (!is.null(names)) {
    cat_df <- cat_df[cat_df$short_name %in% names, , drop = FALSE]
  }
  if (!is.null(domain)) {
    cat_df <- cat_df[cat_df$application_domain %in% domain, , drop = FALSE]
  }
  if (!is.null(bands)) {
    if (inherits(bands, "svid_data")) bands <- bands$bands
    if (inherits(bands, "svid_sensor")) bands <- bands$symbols
    ok <- vapply(strsplit(cat_df$band_list, ","), function(b) all(b %in% bands),
                 logical(1))
    cat_df <- cat_df[ok, , drop = FALSE]
  }
  rownames(cat_df) <- NULL
  if (as_pool) {
    return(.svid_pool(cat_df$formula, NULL, c(catalog = nrow(cat_df)),
                      source = stats::setNames(cat_df$short_name, cat_df$formula)))
  }
  cat_df
}

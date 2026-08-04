#' svid: Systematic Spectral Index Discovery, Deduplication and Validation
#'
#' @description
#' `svid` turns a labelled image time series into a ranked, statistically
#' defensible set of spectral indices. The workflow is deliberately linear:
#' generation, redundancy reduction, fast screening, evaluation, inference,
#' robustness, mapping and reporting.
#'
#' @section Why deduplication comes first:
#' An exhaustively generated pool of arithmetic band combinations contains far
#' fewer distinguishable features than expressions. Two mechanisms are at work.
#' First, different expressions can be algebraically identical: `(N/B)/R` and
#' `(N/R)/B` are both `N/(B*R)`. Second, a rank-based learner such as a random
#' forest cannot distinguish an index from any strictly monotone transformation
#' of it, so `N/R`, `R/N` and `(N-R)/(N+R)` carry identical information for
#' such a model. Ranking a pool without collapsing these classes produces
#' apparent accuracy differences between formulations of one and the same
#' feature and inflates the number of independent candidates being compared.
#' [svid_prune()] performs both reductions and records an audit trail.
#'
#' @section Band symbols:
#' Band symbols follow the Awesome Spectral Indices convention (Montero et al.
#' 2023 \doi{10.1038/s41597-023-02096-0}), so index expressions are portable
#' across sensors and interoperable with the bundled catalogue. See
#' [svid_bands()] and [svid_sensor()].
#'
#' @name svid-package
#' @aliases svid
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  "stage", "n", "label", "index", "value", "metric", "estimate", "lower",
  "upper", "class_name", "epoch", "band", "mean_value", "sd_value", "size",
  "spread", "oa", "tied", "family", "x", "y", "seed", "fold", "score",
  "complexity", "front", "sensor", "delta", "id", "n_members", "canonical",
  "term", "wl", "reflectance", "group", "ymin", "ymax", "count", "a", "b"
))

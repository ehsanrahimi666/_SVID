# ---------------------------------------------------------------------------
# Module 10: graphics
# ---------------------------------------------------------------------------

.gg <- function() {
  .need("ggplot2")
  invisible(TRUE)
}

.svid_colours <- function(n) {
  if (n <= length(.svid_pal)) return(.svid_pal[seq_len(n)])
  grDevices::colorRampPalette(.svid_pal)(n)
}

.svid_theme <- function(base_size = 11) {
  .gg()
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey30", size = base_size - 1),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom")
}

.svid_pal <- c("#2B6A99", "#C0603A", "#4E8A5B", "#8A6BA8", "#B08D3F", "#5E5E5E",
               "#7FA8C9", "#D89A7A", "#8FBF9A", "#B9A5CC", "#D6C084", "#9E9E9E")

#' Redundancy of the generated pool
#'
#' Stage-by-stage reduction from raw expressions to distinguishable features.
#'
#' @param pool A pruned [svid_pool].
#' @param title,subtitle Plot labels.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE)) svid_plot_redundancy(p)
svid_plot_redundancy <- function(pool, title = "Redundancy of the index pool",
                                 subtitle = NULL) {
  .gg()
  df <- svid_redundancy(pool)
  df$stage <- factor(gsub("_", " ", df$stage), levels = gsub("_", " ", df$stage))
  if (is.null(subtitle)) {
    subtitle <- sprintf("%d expressions encode %d distinguishable features (%.0f%% redundant)",
                        df$n[1], df$n[nrow(df)], 100 * (1 - df$n[nrow(df)] / df$n[1]))
  }
  ggplot2::ggplot(df, ggplot2::aes(x = stage, y = n)) +
    ggplot2::geom_col(fill = .svid_pal[1], width = 0.65) +
    ggplot2::geom_text(ggplot2::aes(label = format(n, big.mark = ",")),
                       vjust = -0.4, size = 3.4, fontface = "bold") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL,
                  y = "Number of expressions") +
    .svid_theme()
}

#' Members of monotone-equivalence classes
#'
#' Shows the alternative formulations of one feature side by side. When scores
#' are supplied, differences between bars are by construction artefacts.
#'
#' @param pool A pruned [svid_pool].
#' @param which Name of the class to display; defaults to the largest.
#' @param scores Optional named numeric vector of scores per expression.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE)) svid_plot_equivalence(p)
svid_plot_equivalence <- function(pool, which = NULL, scores = NULL) {
  .gg()
  mem <- attr(pool, "members")
  if (is.null(mem)) stop("this pool has no equivalence information; run svid_prune().",
                         call. = FALSE)
  if (is.null(which)) which <- names(mem)[which.max(lengths(mem))]
  ex <- mem[[which]]
  df <- data.frame(expression = ex,
                   score = if (is.null(scores)) seq_along(ex) else scores[ex],
                   stringsAsFactors = FALSE)
  df$expression <- factor(df$expression, levels = df$expression[order(df$score)])
  ggplot2::ggplot(df, ggplot2::aes(x = score, y = expression)) +
    ggplot2::geom_col(fill = .svid_pal[2], width = 0.6) +
    ggplot2::labs(title = "Equivalent formulations of one feature",
                  subtitle = paste0("canonical representative: ", which,
                                    "  |  ", length(ex), " formulations"),
                  x = if (is.null(scores)) "index within class" else "score",
                  y = NULL) +
    .svid_theme()
}

#' Distribution of scores across the pool
#'
#' @param x A [svid_fit] or [svid_separability()] result.
#' @param metric Metric to plot.
#' @param highlight Expressions to mark with a vertical line.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_ranking(svid_separability(p, d))
svid_plot_ranking <- function(x, metric = NULL, highlight = "(N - R) / (N + R)") {
  .gg()
  if (inherits(x, "svid_fit")) {
    if (is.null(metric)) metric <- if (x$type == "classification") "oa" else "r2"
    df <- data.frame(expression = x$metrics$index, score = x$metrics[[metric]])
    lab <- metric
  } else {
    df <- data.frame(expression = x$expression, score = x$score)
    lab <- paste("separability:", x$metric[1])
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score)) +
    ggplot2::geom_histogram(bins = 30, fill = .svid_pal[1], colour = "white",
                            linewidth = 0.2) +
    ggplot2::labs(title = "Score distribution across the pool",
                  subtitle = sprintf("%d candidates", nrow(df)),
                  x = lab, y = "Number of features") +
    .svid_theme()
  hl <- df[df$expression %in% highlight, , drop = FALSE]
  if (nrow(hl)) {
    p <- p + ggplot2::geom_vline(data = hl,
                                 ggplot2::aes(xintercept = score),
                                 colour = .svid_pal[2], linewidth = 0.8) +
      ggplot2::annotate("text", x = hl$score[1], y = Inf, hjust = 1.05, vjust = 1.6,
                        label = hl$expression[1], colour = .svid_pal[2], size = 3)
  }
  p
}

#' Ranking stability
#'
#' @param x A [svid_stability()] result.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' pool <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                                max_depth = 1)), d)
#' fit <- svid_evaluate(pool[1:6], d, n_seeds = 2, baselines = "none")
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_stability(svid_stability(fit, d))
svid_plot_stability <- function(x) {
  .gg()
  parts <- list(data.frame(source = "seed-to-seed SD", value = x$seed_sd))
  if (!is.null(x$within)) {
    parts[[2]] <- data.frame(source = "within equivalence class",
                             value = x$within$spread)
  }
  df <- do.call(rbind, parts)
  df <- df[is.finite(df$value), , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(x = source, y = value)) +
    ggplot2::geom_boxplot(width = 0.4, outlier.shape = NA, colour = "grey25") +
    ggplot2::geom_jitter(width = 0.08, alpha = 0.45, colour = .svid_pal[1]) +
    ggplot2::geom_hline(yintercept = x$range_top10, linetype = "dashed",
                        colour = .svid_pal[2]) +
    ggplot2::annotate("text", x = 0.6, y = x$range_top10, vjust = -0.6, hjust = 0,
                      size = 3, colour = .svid_pal[2],
                      label = "spread across the top 10") +
    ggplot2::labs(title = "Noise floor of the ranking",
                  subtitle = paste("metric:", x$metric),
                  x = NULL, y = "difference (points)") +
    .svid_theme()
}

#' Per-class accuracy with confidence intervals
#'
#' @param x A [svid_fit] object.
#' @param n Number of top models to show.
#' @param level Confidence level.
#' @return A `ggplot` object.
#' @export
#' @examples
#' \donttest{
#' d <- svid_demo_data()
#' fit <- svid_evaluate(c("(N - R) / (N + R)", "N / (B * R)"), d,
#'                      n_seeds = 2, baselines = "raw_bands")
#' if (requireNamespace("ggplot2", quietly = TRUE)) svid_plot_perclass(fit)
#' }
svid_plot_perclass <- function(x, n = 10L, level = 0.95) {
  .gg()
  if (x$type != "classification") stop("per-class plot needs a classification fit.",
                                       call. = FALSE)
  m <- utils::head(x$metrics, n)
  cols <- grep("^pa_", names(m), value = TRUE)
  nt <- table(x$truth)
  df <- do.call(rbind, lapply(cols, function(cc) {
    cl <- sub("^pa_", "", cc)
    nn <- as.integer(nt[[cl]])
    p <- m[[cc]] / 100
    z <- stats::qnorm(1 - (1 - level) / 2)
    den <- 1 + z^2 / nn; ctr <- p + z^2 / (2 * nn)
    mrg <- z * sqrt(p * (1 - p) / nn + z^2 / (4 * nn^2))
    data.frame(index = m$index, class_name = sprintf("%s (n = %d)", cl, nn),
               estimate = 100 * p, lower = 100 * (ctr - mrg) / den,
               upper = 100 * (ctr + mrg) / den, stringsAsFactors = FALSE)
  }))
  df$index <- factor(df$index, levels = rev(m$index))
  ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = index)) +
    ggplot2::geom_col(fill = .svid_pal[1], alpha = 0.6, width = 0.65) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper),
                            height = 0.22, colour = "grey20") +
    ggplot2::facet_wrap(~ class_name) +
    ggplot2::labs(title = "Per-class producer's accuracy",
                  subtitle = sprintf("%.0f%% Wilson confidence intervals", 100 * level),
                  x = "producer's accuracy (%)", y = NULL) +
    .svid_theme()
}

#' Class or response profiles through the season
#'
#' @param data A [svid_data] object.
#' @param expression Optional index expression; band profiles are drawn if absent.
#' @param ribbon Draw a +/- 1 SD ribbon?
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_profiles(d, "(N - R) / (N + R)")
svid_plot_profiles <- function(data, expression = NULL, ribbon = TRUE) {
  .gg()
  grp <- if (data$type == "classification") as.character(data$response) else
    cut(data$response, 3, labels = c("low", "mid", "high"))
  mats <- if (is.null(expression)) data$values else
    stats::setNames(svid_eval(expression, data), expression)
  df <- do.call(rbind, lapply(names(mats), function(b) {
    m <- mats[[b]]
    do.call(rbind, lapply(unique(grp), function(g) {
      s <- m[grp == g, , drop = FALSE]
      data.frame(band = b, class_name = g, epoch = seq_len(ncol(m)),
                 mean_value = colMeans(s, na.rm = TRUE),
                 sd_value = apply(s, 2, stats::sd, na.rm = TRUE),
                 stringsAsFactors = FALSE)
    }))
  }))
  df$band <- factor(df$band, levels = names(mats))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = epoch, y = mean_value,
                                        colour = class_name, fill = class_name))
  if (ribbon) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_value - sd_value,
                                               ymax = mean_value + sd_value),
                                  alpha = 0.15, colour = NA)
  }
  p + ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ band, scales = "free_y") +
    ggplot2::scale_colour_manual(values = .svid_colours(length(unique(df$class_name)))) +
    ggplot2::scale_fill_manual(values = .svid_colours(length(unique(df$class_name)))) +
    ggplot2::labs(title = if (is.null(expression)) "Band profiles" else "Index profile",
                  x = "epoch", y = "value", colour = NULL, fill = NULL) +
    .svid_theme()
}

#' Two-band separability matrix
#'
#' The classical contour presentation of the hyperspectral index-search
#' literature: separability of the normalised difference of every band pair.
#'
#' @param data A [svid_data] object.
#' @param metric Passed to [svid_separability()].
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' if (requireNamespace("ggplot2", quietly = TRUE)) svid_plot_bandpair(d)
svid_plot_bandpair <- function(data, metric = "auto") {
  .gg()
  b <- data$bands
  gr <- expand.grid(a = b, b = b, stringsAsFactors = FALSE)
  gr <- gr[gr$a != gr$b, ]
  ex <- sprintf("(%s - %s) / (%s + %s)", gr$a, gr$b, gr$a, gr$b)
  sc <- svid_separability(ex, data, metric = metric)
  gr$score <- sc$score[match(ex, sc$expression)]
  gr$a <- factor(gr$a, levels = b); gr$b <- factor(gr$b, levels = b)
  ggplot2::ggplot(gr, ggplot2::aes(x = a, y = b, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(colours = c("#F2F2F2", "#9EBCD4", .svid_pal[1],
                                              "#1B3F5C"), na.value = "white") +
    ggplot2::labs(title = "Normalised-difference separability by band pair",
                  subtitle = "(a - b) / (a + b)", x = "band a", y = "band b",
                  fill = NULL) +
    .svid_theme() + ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Accuracy against complexity, with the Pareto front
#'
#' @param x A [svid_pareto()] result, or an object accepted by it.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_pareto(svid_pareto(svid_separability(p, d)))
svid_plot_pareto <- function(x) {
  .gg()
  if (!inherits(x, "data.frame") || is.null(x$front)) x <- svid_pareto(x)
  ggplot2::ggplot(x, ggplot2::aes(x = complexity, y = score)) +
    ggplot2::geom_point(ggplot2::aes(colour = front), alpha = 0.7, size = 2) +
    ggplot2::geom_step(data = x[x$front, ], colour = .svid_pal[2],
                       direction = "hv", linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "grey70",
                                            `TRUE` = .svid_pal[2]),
                                 labels = c("dominated", "Pareto front")) +
    ggplot2::labs(title = "Accuracy against expression complexity",
                  x = "band occurrences in the expression", y = "score",
                  colour = NULL) +
    .svid_theme()
}

#' Transferability across domains
#'
#' @param x A [svid_transfer()] result.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data(sites = 2)
#' tr <- svid_transfer(c("(N - R) / (N + R)", "N / (B * R)"), d,
#'                     by = "site", n_seeds = 1)
#' if (requireNamespace("ggplot2", quietly = TRUE)) svid_plot_transfer(tr)
svid_plot_transfer <- function(x) {
  .gg()
  df <- data.frame(expression = rep(x$expression, 2),
                   metric = rep(c("within domain", "transferred"),
                                each = nrow(x)),
                   value = c(x$within, x$transferred),
                   stringsAsFactors = FALSE)
  df$expression <- factor(df$expression, levels = rev(x$expression))
  ggplot2::ggplot(df, ggplot2::aes(x = value, y = expression, fill = metric)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7),
                      width = 0.65) +
    ggplot2::scale_fill_manual(values = .svid_pal[c(1, 2)]) +
    ggplot2::labs(title = "Transferability across domains", x = "score",
                  y = NULL, fill = NULL) +
    .svid_theme()
}

#' Separability screening result
#'
#' @param x A [svid_separability()] result.
#' @param n Number of expressions to show.
#' @return A `ggplot` object.
#' @export
#' @examples
#' d <- svid_demo_data()
#' p <- svid_prune(svid_generate(svid_grammar(c("B", "G", "R", "N"),
#'                                             max_depth = 1)), d)
#' if (requireNamespace("ggplot2", quietly = TRUE))
#'   svid_plot_separability(svid_separability(p, d))
svid_plot_separability <- function(x, n = 20L) {
  .gg()
  df <- utils::head(x[!is.na(x$score), , drop = FALSE], n)
  df$expression <- factor(df$expression, levels = rev(df$expression))
  ggplot2::ggplot(df, ggplot2::aes(x = score, y = expression)) +
    ggplot2::geom_col(fill = .svid_pal[3], width = 0.65) +
    ggplot2::labs(title = "Top candidates by separability",
                  subtitle = paste("metric:", df$metric[1]),
                  x = "score", y = NULL) +
    .svid_theme()
}

#' Classified or index raster
#'
#' @param x A `SpatRaster`.
#' @param title,subtitle Labels.
#' @param discrete Treat values as categories?
#' @return A `ggplot` object.
#' @export
#' @examples
#' svid_sensors()
svid_plot_map <- function(x, title = NULL, subtitle = NULL, discrete = NULL) {
  .gg(); .need("terra")
  df <- terra::as.data.frame(x, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "value"
  if (is.null(discrete)) discrete <- is.factor(df$value) || is.character(df$value)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_raster() + ggplot2::coord_equal() +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL,
                  fill = NULL) +
    .svid_theme() +
    ggplot2::theme(axis.text = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())
  if (discrete) p + ggplot2::scale_fill_manual(values = .svid_colours(length(unique(df$value)))) else
    p + ggplot2::scale_fill_viridis_c()
}

#' @export
plot.svid_pool <- function(x, ...) print(svid_plot_redundancy(x, ...))
#' @export
plot.svid_fit <- function(x, ...) print(svid_plot_ranking(x, ...))
#' @export
plot.svid_data <- function(x, ...) print(svid_plot_profiles(x, ...))
#' @export
plot.svid_stability <- function(x, ...) print(svid_plot_stability(x))
#' @export
plot.svid_screen_result <- function(x, ...) print(svid_plot_separability(x, ...))

# ---------------------------------------------------------------------------
# Module 1b: data objects
# ---------------------------------------------------------------------------

#' Create an svid data object from a table
#'
#' The central data structure of the package. A `svid_data` object holds one
#' matrix per spectral band, each with samples in rows and epochs (acquisition
#' dates) in columns, together with the response and per-sample metadata.
#' A single-date study is simply a `svid_data` with one epoch.
#'
#' @param x A data frame in long format: one row per sample and epoch.
#' @param bands Named character vector mapping band symbols to column names of
#'   `x`, e.g. `c(B = "B02", G = "B03", R = "B04", N = "B08")`. Unnamed vectors
#'   are interpreted as columns whose names are already band symbols.
#' @param response Name of the response column. A factor or character column is
#'   treated as classification, a numeric column as regression.
#' @param id Name of the sample identifier column.
#' @param epoch Name of the epoch (date) column. Use `NULL` for single-date data.
#' @param coords Optional length-2 character vector naming the x and y
#'   coordinate columns, required for spatial resampling and mapping.
#' @param group Optional name of a grouping column (parcel, field, plot) used
#'   by grouped resampling.
#' @param site,year Optional column names used by leave-one-site-out and
#'   leave-one-year-out resampling.
#' @param sensor Optional sensor name or [svid_sensor()] object.
#' @param scale Numeric factor applied to all band values, or `NULL` to take
#'   the sensor's scale factor. Use `1` to leave values untouched.
#'
#' @return An object of class `svid_data`.
#' @export
#' @examples
#' data(andong)
#' d <- svid_from_table(andong,
#'                      bands = c(B = "B02", G = "B03", R = "B04", N = "B08"),
#'                      response = "class", id = "point_id", epoch = "date",
#'                      coords = c("x", "y"), group = "parcel",
#'                      sensor = "sentinel2", scale = 1)
#' d
#' dim(d$values$N)
svid_from_table <- function(x, bands, response, id, epoch = NULL,
                            coords = NULL, group = NULL, site = NULL,
                            year = NULL, sensor = NULL, scale = NULL) {
  stopifnot(is.data.frame(x))
  if (is.null(names(bands))) names(bands) <- bands
  miss <- setdiff(c(unname(bands), response, id, epoch, coords, group, site, year),
                  names(x))
  if (length(miss)) stop("column(s) not found in `x`: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  bad <- setdiff(names(bands), .svid_band_table()$symbol)
  if (length(bad)) stop("unknown band symbol(s): ", paste(bad, collapse = ", "),
                        ". See svid_bands().", call. = FALSE)

  if (!is.null(sensor) && !inherits(sensor, "svid_sensor")) sensor <- svid_sensor(sensor)
  if (is.null(scale)) scale <- if (is.null(sensor)) 1 else sensor$scale

  ids <- as.character(x[[id]])
  if (is.null(epoch)) {
    ep <- rep("e1", nrow(x))
  } else {
    ep <- as.character(x[[epoch]])
  }
  uid <- unique(ids); uep <- sort(unique(ep))
  ri <- match(ids, uid); ci <- match(ep, uep)
  if (anyDuplicated(cbind(ri, ci))) {
    stop("duplicated (id, epoch) combinations found; each sample may appear ",
         "only once per epoch.", call. = FALSE)
  }
  n <- length(uid); k <- length(uep)

  values <- lapply(bands, function(cl) {
    m <- matrix(NA_real_, n, k, dimnames = list(uid, uep))
    m[cbind(ri, ci)] <- as.numeric(x[[cl]]) * scale
    m
  })
  names(values) <- names(bands)

  first <- !duplicated(ri)
  meta <- data.frame(id = uid, stringsAsFactors = FALSE)
  ord <- order(ri[first])
  take <- function(col) if (is.null(col)) NULL else x[[col]][first][ord]
  y <- take(response)
  if (is.character(y)) y <- factor(y)
  meta$x <- if (!is.null(coords)) as.numeric(take(coords[1])) else NA_real_
  meta$y <- if (!is.null(coords)) as.numeric(take(coords[2])) else NA_real_
  meta$group <- if (!is.null(group)) as.character(take(group)) else meta$id
  meta$site  <- if (!is.null(site)) as.character(take(site)) else NA_character_
  meta$year  <- if (!is.null(year)) as.character(take(year)) else NA_character_

  .svid_data(values = values, response = y, meta = meta, epochs = uep,
             sensor = sensor, scale = scale)
}

.svid_data <- function(values, response, meta, epochs, sensor = NULL, scale = 1) {
  type <- if (is.factor(response)) "classification" else "regression"
  structure(list(values = values, response = response, meta = meta,
                 epochs = epochs, bands = names(values), sensor = sensor,
                 scale = scale, type = type,
                 n = nrow(meta), n_epochs = length(epochs)),
            class = "svid_data")
}

#' Read a multi-band, multi-date raster cube
#'
#' Wraps a `terra` `SpatRaster` whose layers are ordered date-by-date, four or
#' more bands per date, into an object that [svid_extract()] and [svid_map()]
#' understand.
#'
#' @param x A `SpatRaster`, a file path, or a character vector of file paths
#'   (one per date, each with all bands).
#' @param sensor Sensor name or [svid_sensor()] object.
#' @param bands Character vector of band symbols in layer order within a date.
#'   Defaults to the first `nl / length(dates)` symbols of the sensor.
#' @param dates Character or Date vector of acquisition dates, one per epoch.
#' @param scale Scale factor; defaults to the sensor's.
#'
#' @return An object of class `svid_cube`.
#' @export
#' @examples
#' f <- system.file("extdata", "andong_cube.tif", package = "svid")
#' if (requireNamespace("terra", quietly = TRUE) && nzchar(f)) {
#'   cube <- svid_cube(f, sensor = "sentinel2",
#'                     bands = c("B", "G", "R", "N"),
#'                     dates = svid_cube_dates(f), scale = 1)
#'   cube
#' }
svid_cube <- function(x, sensor, bands = NULL, dates = NULL, scale = NULL) {
  .need("terra")
  if (!inherits(sensor, "svid_sensor")) sensor <- svid_sensor(sensor)
  if (is.character(x)) x <- terra::rast(x)
  if (!inherits(x, "SpatRaster")) stop("`x` must be a SpatRaster or file path.", call. = FALSE)
  nl <- terra::nlyr(x)
  if (is.null(dates)) dates <- "e1"
  dates <- as.character(dates)
  if (is.null(bands)) {
    nb <- nl / length(dates)
    if (nb != round(nb)) {
      stop("cannot infer band count: ", nl, " layers is not a multiple of ",
           length(dates), " dates. Supply `bands`.", call. = FALSE)
    }
    bands <- sensor$symbols[seq_len(nb)]
  }
  if (nl != length(bands) * length(dates)) {
    stop("layer count (", nl, ") does not equal bands (", length(bands),
         ") x dates (", length(dates), ").", call. = FALSE)
  }
  if (is.null(scale)) scale <- sensor$scale
  structure(list(raster = x, sensor = sensor, bands = bands, dates = dates,
                 scale = scale, n_epochs = length(dates)),
            class = "svid_cube")
}

#' Acquisition dates encoded in raster layer names
#'
#' Convenience helper that pulls `YYYY-MM-DD` or `YYYYMMDD` patterns out of the
#' layer names of a raster and returns the unique dates in order.
#'
#' @param x A `SpatRaster` or file path.
#' @return A character vector of unique dates.
#' @export
#' @examples
#' f <- system.file("extdata", "andong_cube.tif", package = "svid")
#' if (requireNamespace("terra", quietly = TRUE) && nzchar(f)) svid_cube_dates(f)
svid_cube_dates <- function(x) {
  .need("terra")
  if (is.character(x)) x <- terra::rast(x)
  nm <- names(x)
  d <- regmatches(nm, regexpr("(19|20)[0-9]{2}[-_]?[0-9]{2}[-_]?[0-9]{2}", nm))
  d <- gsub("[-_]", "", d)
  d <- ifelse(nchar(d) == 8,
              paste(substr(d, 1, 4), substr(d, 5, 6), substr(d, 7, 8), sep = "-"), NA)
  unique(stats::na.omit(d))
}

#' Reference samples
#'
#' Prepare labelled point or polygon samples for extraction. Polygons are
#' converted to their interior pixels; an optional negative buffer removes
#' edge pixels, which are the usual source of mixed-pixel label noise.
#'
#' @param x An `sf` object, a data frame with coordinate columns, or a file
#'   path readable by `sf::st_read()`.
#' @param label Name of the label column.
#' @param coords For data frames, the names of the x and y columns.
#' @param group Optional column identifying parcels or plots. If absent and `x`
#'   is polygonal, one group is assigned per polygon.
#' @param buffer Negative buffer in metres applied to polygons before taking
#'   interior pixels. `0` disables it.
#' @param crs Optional CRS for data-frame input.
#'
#' @return An object of class `svid_samples`.
#' @export
#' @examples
#' pts <- data.frame(x = runif(20, 0, 100), y = runif(20, 0, 100),
#'                   cls = sample(c("a", "b"), 20, TRUE))
#' s <- svid_samples(pts, label = "cls", coords = c("x", "y"))
#' s
svid_samples <- function(x, label, coords = NULL, group = NULL, buffer = 0,
                         crs = NA) {
  if (is.character(x) && length(x) == 1L && file.exists(x)) {
    .need("sf"); x <- sf::st_read(x, quiet = TRUE)
  }
  if (inherits(x, "sf")) {
    .need("sf")
    geom <- sf::st_geometry_type(x, by_geometry = FALSE)
    if (buffer > 0 && as.character(geom) %in% c("POLYGON", "MULTIPOLYGON")) {
      x <- sf::st_buffer(x, -abs(buffer))
      x <- x[!sf::st_is_empty(x), , drop = FALSE]
    }
    if (is.null(group)) {
      grp <- as.character(seq_len(nrow(x)))
    } else {
      grp <- as.character(x[[group]])
    }
    lab <- x[[label]]
    cc <- suppressWarnings(sf::st_coordinates(sf::st_centroid(sf::st_geometry(x))))
    out <- list(geometry = sf::st_geometry(x), label = lab, group = grp,
                x = cc[, 1], y = cc[, 2], crs = sf::st_crs(x),
                polygon = as.character(geom) %in% c("POLYGON", "MULTIPOLYGON"))
  } else {
    stopifnot(is.data.frame(x))
    if (is.null(coords)) stop("`coords` is required for data-frame input.", call. = FALSE)
    lab <- x[[label]]
    grp <- if (is.null(group)) as.character(seq_len(nrow(x))) else as.character(x[[group]])
    out <- list(geometry = NULL, label = lab, group = grp,
                x = as.numeric(x[[coords[1]]]), y = as.numeric(x[[coords[2]]]),
                crs = crs, polygon = FALSE)
  }
  if (is.character(out$label)) out$label <- factor(out$label)
  out$n <- length(out$label)
  structure(out, class = "svid_samples")
}

#' Extract band time series at reference samples
#'
#' @param cube A [svid_cube()].
#' @param samples A [svid_samples()] object.
#' @param fun Aggregation function for polygon samples.
#' @param na_rm Drop samples with any missing band value?
#'
#' @return A [svid_data] object.
#' @export
#' @examples
#' # see vignette("svid-introduction") for a worked raster example
#' data(andong)
#' head(andong, 3)
svid_extract <- function(cube, samples, fun = mean, na_rm = TRUE) {
  .need("terra")
  stopifnot(inherits(cube, "svid_cube"), inherits(samples, "svid_samples"))
  if (!is.null(samples$geometry)) {
    .need("sf")
    v <- terra::vect(samples$geometry)
    ex <- terra::extract(cube$raster, v, fun = if (samples$polygon) fun else NULL,
                         ID = FALSE, na.rm = TRUE)
  } else {
    ex <- terra::extract(cube$raster,
                         cbind(samples$x, samples$y))
  }
  ex <- as.matrix(ex) * cube$scale
  nb <- length(cube$bands); nd <- length(cube$dates)
  values <- lapply(seq_len(nb), function(b) {
    m <- ex[, seq(b, by = nb, length.out = nd), drop = FALSE]
    dimnames(m) <- list(NULL, cube$dates); m
  })
  names(values) <- cube$bands
  meta <- data.frame(id = as.character(seq_len(nrow(ex))),
                     x = samples$x, y = samples$y, group = samples$group,
                     site = NA_character_, year = NA_character_,
                     stringsAsFactors = FALSE)
  keep <- rep(TRUE, nrow(meta))
  if (na_rm) keep <- Reduce(`&`, lapply(values, function(m) stats::complete.cases(m)))
  values <- lapply(values, function(m) m[keep, , drop = FALSE])
  .svid_data(values, samples$label[keep], meta[keep, , drop = FALSE],
             cube$dates, cube$sensor, cube$scale)
}

#' @export
print.svid_data <- function(x, ...) {
  cat("<svid_data>", x$type, "\n")
  cat("  samples:", x$n, " epochs:", x$n_epochs,
      " bands:", paste(x$bands, collapse = " "), "\n")
  if (x$type == "classification") {
    tb <- table(x$response)
    cat("  classes:", paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                            collapse = "  "), "\n")
    if (length(tb) > 1) {
      cat("  imbalance: largest class",
          sprintf("%.1f%%", 100 * max(tb) / sum(tb)),
          "| one sample of the smallest class =",
          sprintf("%.2f pp\n", 100 / min(tb)))
    }
  } else {
    cat("  response:", sprintf("range [%.4g, %.4g] mean %.4g",
                               min(x$response, na.rm = TRUE),
                               max(x$response, na.rm = TRUE),
                               mean(x$response, na.rm = TRUE)), "\n")
  }
  ng <- length(unique(x$meta$group))
  if (ng < x$n) cat("  groups:", ng, "(effective independent units)\n")
  if (!is.null(x$sensor)) cat("  sensor:", x$sensor$name, "\n")
  invisible(x)
}

#' @export
print.svid_samples <- function(x, ...) {
  cat("<svid_samples>", x$n, if (x$polygon) "polygons" else "points", "\n")
  if (is.factor(x$label)) {
    tb <- table(x$label)
    cat("  labels:", paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                           collapse = "  "), "\n")
  }
  cat("  groups:", length(unique(x$group)), "\n")
  invisible(x)
}

#' @export
print.svid_cube <- function(x, ...) {
  cat("<svid_cube>", x$sensor$platform, "\n")
  cat("  bands:", paste(x$bands, collapse = " "),
      " epochs:", x$n_epochs, "\n")
  cat("  dates:", paste(utils::head(x$dates, 3), collapse = ", "),
      if (x$n_epochs > 3) " ..." else "", "\n")
  invisible(x)
}

.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("package '", pkg, "' is required for this function. ",
         "Install it with install.packages(\"", pkg, "\").", call. = FALSE)
  }
  invisible(TRUE)
}

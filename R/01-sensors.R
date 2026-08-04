# ---------------------------------------------------------------------------
# Module 1a: sensors and band symbols
# ---------------------------------------------------------------------------

.svid_band_table <- function() {
  data.frame(
    symbol = c("A", "B", "G1", "G", "Y", "R", "RE1", "RE2", "RE3", "N", "N2",
               "WV", "S1", "S2", "T1", "T2"),
    name = c("Aerosol", "Blue", "Green 1", "Green", "Yellow", "Red",
             "Red Edge 1", "Red Edge 2", "Red Edge 3", "NIR", "NIR 2",
             "Water Vapour", "SWIR 1", "SWIR 2", "Thermal 1", "Thermal 2"),
    wavelength = c(443, 492, 543, 560, 605, 665, 704, 740, 783, 833, 864,
                   945, 1610, 2190, 10900, 12000),
    domain = c("visible", "visible", "visible", "visible", "visible",
               "visible", "red edge", "red edge", "red edge", "nir", "nir",
               "nir", "swir", "swir", "thermal", "thermal"),
    stringsAsFactors = FALSE
  )
}

.svid_sensor_defs <- function() {
  list(
    sentinel2 = list(
      platform = "Sentinel-2 MSI",
      symbols  = c("A", "B", "G", "R", "RE1", "RE2", "RE3", "N", "N2", "WV", "S1", "S2"),
      native   = c(A = "B01", B = "B02", G = "B03", R = "B04", RE1 = "B05",
                   RE2 = "B06", RE3 = "B07", N = "B08", N2 = "B8A", WV = "B09",
                   S1 = "B11", S2 = "B12"),
      wavelength = c(A = 443, B = 492, G = 560, R = 665, RE1 = 704, RE2 = 740,
                     RE3 = 783, N = 833, N2 = 864, WV = 945, S1 = 1614, S2 = 2202),
      resolution = c(A = 60, B = 10, G = 10, R = 10, RE1 = 20, RE2 = 20,
                     RE3 = 20, N = 10, N2 = 20, WV = 60, S1 = 20, S2 = 20),
      scale = 1e-4
    ),
    landsat8 = list(
      platform = "Landsat 8/9 OLI-TIRS",
      symbols  = c("A", "B", "G", "R", "N", "S1", "S2", "T1"),
      native   = c(A = "SR_B1", B = "SR_B2", G = "SR_B3", R = "SR_B4",
                   N = "SR_B5", S1 = "SR_B6", S2 = "SR_B7", T1 = "ST_B10"),
      wavelength = c(A = 443, B = 482, G = 561, R = 655, N = 865,
                     S1 = 1609, S2 = 2201, T1 = 10895),
      resolution = c(A = 30, B = 30, G = 30, R = 30, N = 30,
                     S1 = 30, S2 = 30, T1 = 100),
      scale = 2.75e-5
    ),
    landsat7 = list(
      platform = "Landsat 4-7 TM / ETM+",
      symbols  = c("B", "G", "R", "N", "S1", "S2", "T1"),
      native   = c(B = "SR_B1", G = "SR_B2", R = "SR_B3", N = "SR_B4",
                   S1 = "SR_B5", S2 = "SR_B7", T1 = "ST_B6"),
      wavelength = c(B = 485, G = 560, R = 660, N = 830, S1 = 1650,
                     S2 = 2215, T1 = 11450),
      resolution = c(B = 30, G = 30, R = 30, N = 30, S1 = 30, S2 = 30, T1 = 60),
      scale = 2.75e-5
    ),
    modis = list(
      platform = "MODIS (MOD09)",
      symbols  = c("B", "G", "R", "N", "N2", "S1", "S2"),
      native   = c(B = "sur_refl_b03", G = "sur_refl_b04", R = "sur_refl_b01",
                   N = "sur_refl_b02", N2 = "sur_refl_b05",
                   S1 = "sur_refl_b06", S2 = "sur_refl_b07"),
      wavelength = c(B = 469, G = 555, R = 645, N = 858, N2 = 1240,
                     S1 = 1640, S2 = 2130),
      resolution = c(B = 500, G = 500, R = 250, N = 250, N2 = 500,
                     S1 = 500, S2 = 500),
      scale = 1e-4
    ),
    planetscope = list(
      platform = "PlanetScope SuperDove",
      symbols  = c("A", "B", "G1", "G", "Y", "R", "RE1", "N"),
      native   = c(A = "coastal_blue", B = "blue", G1 = "green_i", G = "green",
                   Y = "yellow", R = "red", RE1 = "rededge", N = "nir"),
      wavelength = c(A = 443, B = 490, G1 = 531, G = 565, Y = 610, R = 665,
                     RE1 = 705, N = 865),
      resolution = c(A = 3, B = 3, G1 = 3, G = 3, Y = 3, R = 3, RE1 = 3, N = 3),
      scale = 1e-4
    ),
    worldview3 = list(
      platform = "WorldView-3",
      symbols  = c("A", "B", "G", "Y", "R", "RE1", "N", "N2", "S1", "S2"),
      native   = c(A = "coastal", B = "blue", G = "green", Y = "yellow",
                   R = "red", RE1 = "rededge", N = "nir1", N2 = "nir2",
                   S1 = "swir3", S2 = "swir6"),
      wavelength = c(A = 425, B = 480, G = 545, Y = 605, R = 660, RE1 = 725,
                     N = 833, N2 = 950, S1 = 1660, S2 = 2205),
      resolution = c(A = 1.24, B = 1.24, G = 1.24, Y = 1.24, R = 1.24,
                     RE1 = 1.24, N = 1.24, N2 = 1.24, S1 = 3.7, S2 = 3.7),
      scale = 1
    )
  )
}

#' Standard spectral band symbols
#'
#' Returns the band-symbol vocabulary used throughout `svid`. Symbols follow
#' the Awesome Spectral Indices convention so that index expressions are
#' portable between sensors.
#'
#' @param sensor Optional sensor name (see [svid_sensors()]). If supplied,
#'   only the symbols available on that sensor are returned.
#' @param domain Optional character vector restricting the result to one or
#'   more spectral domains: `"visible"`, `"red edge"`, `"nir"`, `"swir"`,
#'   `"thermal"`.
#'
#' @return A data frame with columns `symbol`, `name`, `wavelength` (nm) and
#'   `domain`.
#' @export
#' @examples
#' svid_bands()
#' svid_bands(sensor = "landsat8")
#' svid_bands(domain = c("nir", "swir"))
svid_bands <- function(sensor = NULL, domain = NULL) {
  tb <- .svid_band_table()
  if (!is.null(sensor)) {
    s <- svid_sensor(sensor)
    tb <- tb[tb$symbol %in% s$symbols, , drop = FALSE]
    tb$wavelength <- unname(s$wavelength[tb$symbol])
  }
  if (!is.null(domain)) {
    domain <- match.arg(domain, unique(tb$domain), several.ok = TRUE)
    tb <- tb[tb$domain %in% domain, , drop = FALSE]
  }
  rownames(tb) <- NULL
  tb
}

#' Names of the built-in sensor definitions
#'
#' @return A character vector of sensor keys accepted by [svid_sensor()].
#' @export
#' @examples
#' svid_sensors()
svid_sensors <- function() names(.svid_sensor_defs())

#' Sensor definition
#'
#' Look up a built-in sensor, or construct a custom one. A sensor definition
#' maps standard band symbols onto the native band names of a platform and
#' records central wavelengths, spatial resolutions and the scale factor that
#' converts stored digital numbers to reflectance.
#'
#' @param name Sensor key (see [svid_sensors()]), or a name for the custom
#'   sensor being defined.
#' @param symbols Character vector of band symbols, for custom sensors.
#' @param native Named character vector mapping symbols to native band names.
#' @param wavelength Named numeric vector of central wavelengths in nm.
#' @param resolution Named numeric vector of pixel sizes in metres.
#' @param scale Numeric scale factor converting stored values to reflectance
#'   in \[0, 1\]. Use `1` if the data are already reflectance.
#' @param platform Human-readable platform description.
#'
#' @return An object of class `svid_sensor`.
#' @export
#' @examples
#' s2 <- svid_sensor("sentinel2")
#' s2
#' s2$native[["S1"]]
#'
#' # a custom four-band UAV camera
#' uav <- svid_sensor("micasense", symbols = c("B", "G", "R", "RE1", "N"),
#'                    native = c(B = "b1", G = "b2", R = "b3",
#'                               RE1 = "b4", N = "b5"),
#'                    wavelength = c(B = 475, G = 560, R = 668,
#'                                   RE1 = 717, N = 840),
#'                    resolution = rep(0.05, 5), scale = 1)
#' uav
svid_sensor <- function(name, symbols = NULL, native = NULL, wavelength = NULL,
                        resolution = NULL, scale = 1, platform = NULL) {
  defs <- .svid_sensor_defs()
  if (is.null(symbols) && name %in% names(defs)) {
    out <- defs[[name]]
    out$name <- name
  } else {
    if (is.null(symbols)) {
      stop("Unknown sensor '", name, "'. Available: ",
           paste(names(defs), collapse = ", "),
           ". Supply `symbols` to define a custom sensor.", call. = FALSE)
    }
    bad <- setdiff(symbols, .svid_band_table()$symbol)
    if (length(bad)) {
      stop("Unknown band symbol(s): ", paste(bad, collapse = ", "),
           ". See svid_bands() for the accepted vocabulary.", call. = FALSE)
    }
    if (is.null(native)) native <- stats::setNames(symbols, symbols)
    if (is.null(wavelength)) {
      tb <- .svid_band_table()
      wavelength <- stats::setNames(tb$wavelength[match(symbols, tb$symbol)], symbols)
    }
    if (is.null(resolution)) resolution <- stats::setNames(rep(NA_real_, length(symbols)), symbols)
    out <- list(name = name, platform = platform %||% name, symbols = symbols,
                native = native, wavelength = wavelength,
                resolution = resolution, scale = scale)
  }
  structure(out, class = "svid_sensor")
}

#' @export
print.svid_sensor <- function(x, ...) {
  cat("<svid_sensor>", x$name, "--", x$platform, "\n")
  cat("  bands (", length(x$symbols), "): ",
      paste(x$symbols, collapse = " "), "\n", sep = "")
  cat("  scale factor: ", format(x$scale), "\n", sep = "")
  tb <- data.frame(symbol = x$symbols,
                   native = unname(x$native[x$symbols]),
                   nm = unname(x$wavelength[x$symbols]),
                   m = unname(x$resolution[x$symbols]))
  print(tb, row.names = FALSE)
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

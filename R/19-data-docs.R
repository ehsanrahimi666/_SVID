#' Sentinel-2A time series for three land covers, Andong, South Korea
#'
#' Surface-reflectance time series extracted at 1,790 field-verified reference
#' points in the Angyo agricultural area of Andong City, North Gyeongsang
#' Province, Republic of Korea. The site covers 323.9 ha between 457,320 and
#' 459,690 m easting and 4,045,650 and 4,048,120 m northing in UTM zone 52N
#' (EPSG:32652), equivalently 36.556-36.600 N and 128.522-128.573 E, and falls
#' entirely within Sentinel-2 tile T52SDF.
#'
#' @format A data frame with 57,280 rows (1,790 points x 32 dates) and 11
#'   columns:
#' \describe{
#'   \item{point_id}{Reference point identifier.}
#'   \item{split}{`training` (1,290 points digitised as contiguous blocks
#'     inside 49 field parcels) or `validation` (500 independently located
#'     points withheld from model fitting).}
#'   \item{class}{Land cover: `greenhouse`, `pepper` or `rice`. Training counts
#'     are 212, 173 and 905; validation counts are 58, 13 and 429.}
#'   \item{x, y}{Easting and northing in metres, EPSG:32652.}
#'   \item{date}{Acquisition date.}
#'   \item{B02, B03, B04, B08}{Blue, green, red and near-infrared surface
#'     reflectance as Level-2A digital numbers at 10 m; divide by 10,000 for
#'     reflectance.}
#'   \item{parcel}{Contiguous digitised block, obtained by single-linkage
#'     clustering of the coordinates at a 15 m threshold. The 1,290 training
#'     pixels fall into 49 blocks; the 500 validation points are spatially
#'     isolated and each forms its own block, giving 491 levels in total.
#'     Grouped resampling should use this column rather than treating pixels as
#'     independent.}
#' }
#'
#' @details
#' Imagery are Sentinel-2A Level-2A products for the 2022 calendar year,
#' atmospherically corrected with Sen2Cor, cloud- and shadow-masked using the
#' QA60 band, and restricted to the four 10 m bands. Thirty-two acquisitions
#' were cloud-free over the site: three in January, three in February, three in
#' March, four in April, five in May, two in June, one in July, two in August,
#' three in October, five in November and one in December. No cloud-free
#' acquisition was available in September.
#'
#' Reference labels were assigned by field survey where parcels were
#' accessible and by interpretation of high-resolution imagery elsewhere. The
#' three classes are phenologically distinct: rice paddies fall to an NDVI near
#' zero in mid-May when fields are flooded before transplanting and peak at
#' about 0.61 on 22 August; pepper remains below 0.21 until June and peaks at
#' about 0.69 in late August; greenhouses show no seasonal cycle, remaining
#' between 0.10 and 0.33 all year because the plastic cover dominates the
#' signal.
#'
#' The validation set is strongly imbalanced, rice accounting for 85.8\% and
#' pepper for 2.6\% of the 500 points. One pepper sample is worth 7.7 percentage
#' points of that class's producer's accuracy. The data are bundled partly for
#' this reason: they are a realistic rather than a convenient benchmark.
#'
#' @source Sentinel-2 Level-2A imagery, Copernicus Open Access Hub / Copernicus
#'   Data Space Ecosystem, tile T52SDF, 2022. Contains modified Copernicus
#'   Sentinel data 2022. Reference points collected by the authors.
#' @examples
#' data(andong)
#' str(andong)
#' table(andong$split[!duplicated(andong$point_id)],
#'       andong$class[!duplicated(andong$point_id)])
#'
#' d <- svid_from_table(andong,
#'                      bands = c(B = "B02", G = "B03", R = "B04", N = "B08"),
#'                      response = "class", id = "point_id", epoch = "date",
#'                      coords = c("x", "y"), group = "parcel",
#'                      sensor = "sentinel2")
#' d
"andong"

#' Catalogue of published spectral indices
#'
#' The Awesome Spectral Indices catalogue: 280 published spectral indices with
#' standardised formulae expressed over the band symbols used throughout this
#' package. Provided so that discovered indices can be benchmarked against the
#' published literature rather than against NDVI alone.
#'
#' @format A data frame with 280 rows and 8 columns:
#' \describe{
#'   \item{short_name}{Conventional acronym, e.g. `NDVI`.}
#'   \item{long_name}{Full name.}
#'   \item{formula}{Formula over standard band symbols, directly evaluable by
#'     [svid_eval()].}
#'   \item{band_list}{Comma-separated band symbols required.}
#'   \item{application_domain}{`vegetation`, `water`, `burn`, `urban`, `soil`,
#'     `snow`, `radar`, `clouds` or `kernel`.}
#'   \item{platforms}{Semicolon-separated sensors on which the index is defined.}
#'   \item{reference}{DOI or URL of the original publication.}
#'   \item{date_of_addition}{Date the entry was added to the catalogue.}
#' }
#' @source Montero, D., Aybar, C., Mahecha, M.D., Martinuzzi, F., Soechting, M.
#'   and Wieneke, S. (2023) A standardized catalogue of spectral indices to
#'   advance the use of remote sensing in Earth system research.
#'   \emph{Scientific Data} 10, 197. \doi{10.1038/s41597-023-02096-0}.
#'   Catalogue retrieved from
#'   \url{https://github.com/awesome-spectral-indices/awesome-spectral-indices}.
#' @examples
#' data(asi_indices)
#' head(asi_indices[, c("short_name", "formula", "application_domain")])
#' subset(asi_indices, short_name %in% c("NDVI", "EVI", "NDWI"))$formula
"asi_indices"

#' Landsat-8 crop-year time series for the Brazilian Cerrado
#'
#' Surface-reflectance time series at 5,163 labelled points across the Brazilian
#' Cerrado biome, covering one complete crop year at sixteen-day cadence. The
#' data support phenology-based class discrimination: the eight classes differ
#' less in their spectral magnitude than in the timing of their green-up and
#' senescence, so a single acquisition separates them poorly while the annual
#' trajectory separates them well.
#'
#' @format A data frame with 123,912 rows (5,163 points x 24 dates) and 12
#'   columns:
#' \describe{
#'   \item{point_id}{Reference point identifier.}
#'   \item{class}{Land cover, one of `Annual_Crop`, `Cerradao`, `Cerrado`,
#'     `Open_Cerrado`, `Pasture`, `Silviculture`, `Sugarcane` (700 points each)
#'     and `Water` (263 points).}
#'   \item{longitude, latitude}{Decimal degrees, WGS84. The samples span
#'     41.7-60.1 W and 2.4-24.6 S.}
#'   \item{date}{Composite date, 24 sixteen-day steps from 2017-08-29 to
#'     2018-08-29.}
#'   \item{SR_B1 ... SR_B7}{Landsat-8 OLI surface reflectance in the coastal,
#'     blue, green, red, near-infrared, SWIR-1 and SWIR-2 bands, already scaled
#'     to the 0-1 reflectance range.}
#'   \item{tile}{One-degree tile identifier, used as the grouping variable for
#'     spatially blocked resampling.}
#' }
#'
#' @details
#' The source samples were collected for the Brazil Data Cube and distributed
#' with the `sitsdata` package as `samples_cerrado_lc8` (48,850 points). The
#' version bundled here drops classes with fewer than 200 points and takes a
#' stratified random subsample of at most 700 points per class, so that the
#' data set is small enough to distribute while retaining the full temporal
#' resolution and the natural imbalance of the `Water` class.
#'
#' Band identity was verified rather than assumed: NDVI computed from `SR_B5`
#' and `SR_B4` reproduces the NDVI supplied with the source data to three
#' decimal places, and the pooled mean spectrum shows the expected vegetation
#' shape, with a red minimum near 0.08, a near-infrared plateau near 0.29 and a
#' shortwave-infrared decline from 0.25 to 0.15.
#'
#' @source Brazil Data Cube, Landsat-8 OLI surface reflectance, 2017-2018;
#'   sample set distributed with the `sitsdata` package
#'   (\url{https://github.com/e-sensing/sitsdata}). Described in Picoli, M.C.A.
#'   et al. (2018) Big Earth observation time series analysis for monitoring
#'   Brazilian agriculture, \emph{ISPRS Journal of Photogrammetry and Remote
#'   Sensing} 145, 328-339, and Simoes, R. et al. (2021) Satellite Image Time
#'   Series Analysis for Big Earth Observation Data, \emph{Remote Sensing} 13,
#'   2428. \doi{10.3390/rs13132428}
#' @examples
#' data(cerrado)
#' str(cerrado, max.level = 1)
#' table(cerrado$class[!duplicated(cerrado$point_id)])
#'
#' d <- svid_from_table(cerrado,
#'                      bands = c(A = "SR_B1", B = "SR_B2", G = "SR_B3",
#'                                R = "SR_B4", N = "SR_B5", S1 = "SR_B6",
#'                                S2 = "SR_B7"),
#'                      response = "class", id = "point_id", epoch = "date",
#'                      coords = c("longitude", "latitude"), group = "tile",
#'                      sensor = "landsat8", scale = 1)
#' d
"cerrado"

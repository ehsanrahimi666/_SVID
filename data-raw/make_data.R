## Build the bundled data sets. Run from the package root.
U <- "/mnt/user-data/uploads/"

## ---- andong: Sentinel-2A time series, Andong, South Korea -----------------
p <- read.csv(paste0(U, "andong_reflectance_points.csv"), stringsAsFactors = FALSE)
p <- p[, c("point_id", "split", "class_name", "x", "y", "date",
           "B02", "B03", "B04", "B08")]
names(p)[names(p) == "class_name"] <- "class"
p$class <- factor(p$class, levels = c("greenhouse", "pepper", "rice"))
p$split <- factor(p$split, levels = c("training", "validation"))

## parcel identifier: single-linkage clustering of the coordinates at 15 m,
## which recovers the contiguous digitised blocks
first <- !duplicated(p$point_id)
xy <- as.matrix(p[first, c("x", "y")])
d <- dist(xy)
h <- hclust(d, method = "single")
parcel <- cutree(h, h = 15)
p$parcel <- sprintf("P%03d", parcel[match(p$point_id, p$point_id[first])])

andong <- p[order(p$point_id, p$date), ]
rownames(andong) <- NULL
usethis_save <- function(obj, name) {
  assign(name, obj)
  save(list = name, file = file.path("data", paste0(name, ".rda")),
       compress = "xz", version = 2)
}
usethis_save(andong, "andong")

## ---- asi_indices: Awesome Spectral Indices catalogue ----------------------
j <- jsonlite_read <- function(f) {
  txt <- paste(readLines(f, warn = FALSE), collapse = "")
  ## minimal JSON reader is avoided: use the python-produced flat CSV instead
  txt
}
asi <- read.csv("data-raw/asi_flat.csv", stringsAsFactors = FALSE)
asi_indices <- asi[order(asi$short_name), ]
rownames(asi_indices) <- NULL
usethis_save(asi_indices, "asi_indices")

cat("andong:", nrow(andong), "rows;  asi_indices:", nrow(asi_indices), "rows\n")

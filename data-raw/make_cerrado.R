## Build the bundled `cerrado` data set from sitsdata::samples_cerrado_lc8
## Run from the package root; set CERRADO_RDA if samples_cerrado_lc8.rda is elsewhere.
e <- new.env(); load(Sys.getenv("CERRADO_RDA", file.path("data-raw", "samples_cerrado_lc8.rda")), envir = e)
x <- get(ls(e)[1], e)

## drop classes with too few samples to support an interval estimate
tb <- table(x$label)
x <- x[x$label %in% names(tb)[tb >= 200], ]
x$label <- factor(as.character(x$label))

## stratified subsample: at most 700 points per class, keeps the package small
set.seed(42)
idx <- unlist(lapply(split(seq_len(nrow(x)), x$label),
                     function(i) if (length(i) > 700) sample(i, 700) else i))
x <- x[sort(idx), ]

bands <- c("B1", "B2", "B3", "B4", "B5", "B6", "B7")
rows <- lapply(seq_len(nrow(x)), function(i) {
  ts <- x$time_series[[i]]
  data.frame(point_id = sprintf("c%05d", i),
             class = as.character(x$label[i]),
             longitude = x$longitude[i], latitude = x$latitude[i],
             date = as.character(ts$Index),
             ts[, bands], stringsAsFactors = FALSE)
})
cerrado <- do.call(rbind, rows)
names(cerrado)[names(cerrado) %in% bands] <-
  c("SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6", "SR_B7")
cerrado$class <- factor(cerrado$class)
## spatial blocks for grouped resampling: one-degree tiles
cerrado$tile <- sprintf("T%+04d%+04d", floor(cerrado$longitude),
                        floor(cerrado$latitude))
rownames(cerrado) <- NULL
save(cerrado, file = "data/cerrado.rda", compress = "xz", version = 2)
cat("cerrado:", nrow(cerrado), "rows |", length(unique(cerrado$point_id)), "points |",
    length(unique(cerrado$date)), "dates |", nlevels(cerrado$class), "classes |",
    length(unique(cerrado$tile)), "tiles\n")
print(table(cerrado$class[!duplicated(cerrado$point_id)]))

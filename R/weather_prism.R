PRISM_URL <- "https://services.nacse.org/prism/data/get/us/4km/%s/%s"
PRISM_VARS <- c(ppt = "RAIN", tmax = "TMAX", tmin = "TMIN", tdmean = "TDEW")
# Polite spacing between NACSE requests (seconds) to avoid throttle responses.
PRISM_REQUEST_DELAY <- 1.0

.prism_is_zip <- function(path) {
  # A valid zip starts with the magic bytes "PK" (0x50 0x4B). NACSE throttle
  # responses are small HTML/text bodies that fail this check.
  if (!file.exists(path) || file.info(path)$size < 4) return(FALSE)
  magic <- readBin(path, "raw", n = 2)
  identical(magic, as.raw(c(0x50, 0x4B)))
}

.prism_download_grid <- function(var, day, cache_dir) {
  ymd <- format(day, "%Y%m%d")
  out_dir <- file.path(cache_dir, var, ymd)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  existing <- list.files(out_dir, pattern = "\\.(bil|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  if (length(existing)) return(existing[1])
  zip_path <- file.path(out_dir, paste0(var, "_", ymd, ".zip"))
  ok <- tryCatch({
    Sys.sleep(PRISM_REQUEST_DELAY)
    utils::download.file(sprintf(PRISM_URL, var, ymd), zip_path, mode = "wb", quiet = TRUE)
    if (.prism_is_zip(zip_path)) {
      utils::unzip(zip_path, exdir = out_dir)
      TRUE
    } else {
      message(sprintf("  PRISM %s %s: response was not a valid zip (likely throttled); skipping.", var, ymd))
      FALSE
    }
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) return(NA_character_)
  existing <- list.files(out_dir, pattern = "\\.(bil|tif|tiff)$", full.names = TRUE, ignore.case = TRUE)
  if (length(existing)) existing[1] else NA_character_
}

process_weather_prism <- function(shapefile, start_year, end_year, output_dir,
                                  id_col, lat_col, lon_col, n_cores, log_file,
                                  prism_cache_dir) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("package 'terra' required for PRISM")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(prism_cache_dir, recursive = TRUE, showWarnings = FALSE)
  dates <- seq(as.Date(sprintf("%d-01-01", start_year)),
               min(as.Date(sprintf("%d-12-31", end_year)), Sys.Date() - 2),
               by = "day")
  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)
  frames <- setNames(vector("list", length(ids)), ids)
  for (i in seq_along(frames)) frames[[i]] <- list()
  message(sprintf("--- Starting PRISM Processing (Years: %d-%d) ---", start_year, end_year))
  # NB: iterate by index. `for (day in dates)` coerces each Date element to its
  # underlying integer, which breaks format(day, ...); dates[di] preserves class.
  for (di in seq_along(dates)) {
    day <- dates[di]
    day_vals <- list()
    for (var in names(PRISM_VARS)) {
      p <- .prism_download_grid(var, day, prism_cache_dir)
      if (!is.na(p)) day_vals[[PRISM_VARS[[var]]]] <- as.numeric(terra::extract(terra::rast(p), pts_vect, ID = FALSE)[, 1])
    }
    for (i in seq_along(ids)) {
      frames[[i]][[length(frames[[i]]) + 1]] <- data.frame(
        DATE = sprintf("%d%03d", as.integer(format(day, "%Y")), as.integer(format(day, "%j"))),
        YEAR = as.integer(format(day, "%Y")), MM = as.integer(format(day, "%m")),
        SRAD = -99, TMAX = if (!is.null(day_vals$TMAX)) day_vals$TMAX[i] else NA,
        TMIN = if (!is.null(day_vals$TMIN)) day_vals$TMIN[i] else NA,
        RAIN = if (!is.null(day_vals$RAIN)) day_vals$RAIN[i] else NA,
        TDEW = if (!is.null(day_vals$TDEW)) day_vals$TDEW[i] else -99,
        RH2M = -99, WIND = -99)
    }
  }
  written <- 0
  for (i in seq_along(ids)) {
    df <- do.call(rbind, frames[[i]])
    df <- df[!is.na(df$TMAX) & !is.na(df$TMIN), ]
    if (!nrow(df)) {
      write(sprintf("PRISM point %s: no valid TMAX/TMIN data extracted", ids[i]), file = log_file, append = TRUE)
      next
    }
    weather_write_wth(df, ids[i], lats[i], lons[i], output_dir, "PRISM 4km", "PRSM", wndht = -99)
    written <- written + 1
  }
  message(sprintf("\nPRISM processing complete: %d/%d point(s) written.\n", written, length(ids)))
}

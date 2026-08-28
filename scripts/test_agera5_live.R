#!/usr/bin/env Rscript

# Live AgERA5 smoke test for the R implementation.
#
# Uses the exact point at the center of the bounding box seen in the Bioenergy
# download log: 33.6816 N, -102.5220 E.  The default test downloads only one
# AgERA5 variable/year (TMIN, 1987), validates the canonical ZIP cache, opens the
# returned NetCDF with terra, and extracts the real value at that point.
#
# Requirements:
#   - a valid Copernicus CDS Personal Access Token configured for ecmwfr /
#     dssatutils (same credentials used by the production pipeline)
#   - packages ecmwfr and terra
#
# Run from the dssatutils repository:
#   Rscript scripts/test_agera5_live.R

suppressPackageStartupMessages({
  library(dssatutils)
  library(ecmwfr)
  library(terra)
})

lat <- 33.6816
lon <- -102.5220
year <- 1987L
cache_dir <- Sys.getenv(
  "AGERA5_LIVE_TEST_CACHE",
  unset = file.path(tempdir(), "dssatutils_agera5_live")
)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("dssatutils version: %s\n", as.character(utils::packageVersion("dssatutils"))))
pd <- utils::packageDescription("dssatutils")
remote_sha <- if (!is.null(pd$RemoteSha)) pd$RemoteSha else "<not recorded>"
cat(sprintf("dssatutils RemoteSha: %s\n", remote_sha))
cat(sprintf("Live point: lat=%.4f lon=%.4f; year=%d\n", lat, lon, year))
cat(sprintf("Cache: %s\n", cache_dir))

# Fail early with a useful credentials error rather than after building a job.
dssatutils:::.agera5_ensure_ecmwfr_key(quiet = FALSE)

area <- c(lat + 0.2, lon - 0.2, lat - 0.2, lon + 0.2) # N,W,S,E
spec <- dssatutils:::.agera5_vars$TMIN
job <- dssatutils:::.agera5_job("TMIN", year, spec, area, cache_dir)

cat(sprintf("Expected canonical archive: %s\n", job$zip_dest))
cat("Submitting one real AgERA5 variable-year request...\n")
res <- dssatutils:::.agera5_download_job(job)
cat(res$message, "\n", sep = "")

if (!isTRUE(res$ok)) stop(res$message, call. = FALSE)
if (!file.exists(job$zip_dest)) {
  stop("Download reported success but canonical .zip cache does not exist: ", job$zip_dest,
       call. = FALSE)
}
legacy_double <- paste0(job$zip_dest, ".zip")
if (file.exists(legacy_double)) {
  stop("Legacy .zip.zip artifact still exists after download/recovery: ", legacy_double,
       call. = FALSE)
}
if (!dssatutils:::.agera5_valid_zip(job$zip_dest)) {
  stop("Canonical AgERA5 archive failed ZIP/NetCDF validation: ", job$zip_dest,
       call. = FALSE)
}

nc_files <- res$data_files
if (!length(nc_files)) stop("No NetCDF file was exposed from the archive.", call. = FALSE)
cat(sprintf("NetCDF file(s): %s\n", paste(basename(nc_files), collapse = ", ")))

r <- suppressWarnings(terra::rast(nc_files))
if (terra::nlyr(r) < 1L) stop("AgERA5 NetCDF contains no raster layers.", call. = FALSE)
pts <- terra::vect(data.frame(lon = lon, lat = lat), geom = c("lon", "lat"), crs = "EPSG:4326")
ex <- suppressWarnings(terra::extract(r, pts, ID = FALSE))
vals_k <- as.numeric(ex[1, ])
finite <- vals_k[is.finite(vals_k)]
if (!length(finite)) stop("No finite AgERA5 TMIN value extracted at the live point.", call. = FALSE)

vals_c <- finite - 273.15
cat(sprintf("Extracted %d finite daily TMIN value(s).\n", length(finite)))
cat(sprintf("First finite TMIN: %.2f K = %.2f C\n", finite[1], vals_c[1]))
cat(sprintf("TMIN range: %.2f to %.2f C\n", min(vals_c), max(vals_c)))
cat("PASS: real AgERA5 R download, canonical cache, NetCDF open, and point extraction all succeeded.\n")

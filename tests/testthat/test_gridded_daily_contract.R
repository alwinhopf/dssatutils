test_that("NetCDF daily adapter rejects duplicate and subdaily records", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = 2, xmin = 9, xmax = 11,
                   ymin = 39, ymax = 41, crs = "EPSG:4326")
  terra::values(r) <- matrix(c(10, 10), nrow = 1)
  pts <- terra::vect(data.frame(x = 10, y = 40), geom = c("x", "y"), crs = "EPSG:4326")
  path <- tempfile(fileext = ".tif")
  terra::time(r) <- as.POSIXct(c("2001-01-01 00:00:00", "2001-01-01 00:00:00"), tz = "UTC")
  terra::writeRaster(r, path, overwrite = TRUE)
  expect_error(weather_extract_netcdf_series(path, "p", pts, 2001, 2001, "rain"),
               "Duplicate weather timestamps")
  terra::time(r) <- as.POSIXct(c("2001-01-01 00:00:00", "2001-01-01 12:00:00"), tz = "UTC")
  terra::writeRaster(r, path, overwrite = TRUE)
  expect_error(weather_extract_netcdf_series(path, "p", pts, 2001, 2001, "rain"),
               "Subdaily weather data")
  terra::time(r) <- as.Date(c("2001-01-01", "2001-01-02"))
  terra::values(r) <- matrix(c(1, 2) / 86400, nrow = 1)
  terra::units(r) <- "kg m-2 s-1"
  terra::writeRaster(r, path, overwrite = TRUE, datatype = "FLT8S")
  result <- weather_extract_netcdf_series(path, "p", pts, 2001, 2001, "rain")
  expect_equal(result$p, c("2001001" = 1, "2001002" = 2))
  pts_outside <- terra::vect(data.frame(x = 50, y = 80), geom = c("x", "y"), crs = "EPSG:4326")
  expect_error(weather_extract_netcdf_series(path, "p_out", pts_outside, 2001, 2001, "rain"),
               "outside grid domain")
})

test_that("NetCDF extraction handles 0..360 and -180..180 longitude coordinate wrapping", {
  skip_if_not_installed("terra")
  path <- tempfile(fileext = ".tif")

  # 1. Raster in 0..360: xmin=279, xmax=281, query point at x=-80
  r0_360 <- terra::rast(nrows = 1, ncols = 1, nlyrs = 1, xmin = 279, xmax = 281,
                        ymin = 24, ymax = 26, crs = "EPSG:4326")
  terra::time(r0_360) <- as.Date("2001-01-01")
  terra::units(r0_360) <- "mm"
  terra::values(r0_360) <- matrix(5.0, nrow = 1)
  terra::writeRaster(r0_360, path, overwrite = TRUE)
  pts_neg <- terra::vect(data.frame(x = -80, y = 25), geom = c("x", "y"), crs = "EPSG:4326")
  res1 <- weather_extract_netcdf_series(path, "p1", pts_neg, 2001, 2001, "rain")
  expect_equal(res1$p1, c("2001001" = 5.0))

  # 2. Raster in -180..180: xmin=-81, xmax=-79, query point at x=280
  r_180 <- terra::rast(nrows = 1, ncols = 1, nlyrs = 1, xmin = -81, xmax = -79,
                       ymin = 24, ymax = 26, crs = "EPSG:4326")
  terra::time(r_180) <- as.Date("2001-01-01")
  terra::units(r_180) <- "mm"
  terra::values(r_180) <- matrix(7.0, nrow = 1)
  terra::writeRaster(r_180, path, overwrite = TRUE)
  pts_pos <- terra::vect(data.frame(x = 280, y = 25), geom = c("x", "y"), crs = "EPSG:4326")
  res2 <- weather_extract_netcdf_series(path, "p2", pts_pos, 2001, 2001, "rain")
  expect_equal(res2$p2, c("2001001" = 7.0))
})

test_that("weather_convert_units parity and edge cases", {
  # Unknown units fail rather than guessing from magnitude
  expect_error(weather_convert_units(25.0, "unknown", "temp"), "units")

  # Temp: 'degK', 'kelvin', 'K' should trigger subtraction
  expect_equal(weather_convert_units(298.15, "degk", "temp"), 25.0, tolerance = 1e-4)
  expect_equal(weather_convert_units(298.15, "kelvin", "temp"), 25.0, tolerance = 1e-4)

  # Missing radiation units fail explicitly
  expect_error(weather_convert_units(15000000, "", "srad"), "units")

  # SRAD: 'MJ' units preserved
  expect_equal(weather_convert_units(18.5, "MJ m-2 day-1", "srad"), 18.5)

  # Wind: 10 m -> 2 m
  expect_equal(weather_convert_units(4.0, "m s-1", "wind", wind_height_m = 10), 4.0 * 0.748, tolerance = 1e-3)

  # Vapour pressure to dewpoint
  td <- weather_convert_units(12.27, "hPa", "vp")
  expect_true(abs(td - 10.0) < 0.5)
})


test_that("GridMET AMP is the full mean annual temperature range", {
  dates <- seq(as.Date("2001-01-01"), as.Date("2002-12-31"), by = "day")
  tmean <- as.numeric(format(dates, "%m"))
  wth <- data.frame(DATE = dates, TMAX = tmean + 5, TMIN = tmean - 5)
  expect_equal(dssatutils:::.gridmet_calc_amp(wth), 11)
  tmean[format(dates, "%Y") == "2002"] <- 13 - tmean[format(dates, "%Y") == "2002"]
  wth$TMAX <- tmean + 5
  wth$TMIN <- tmean - 5
  expect_equal(dssatutils:::.gridmet_calc_amp(wth), 11)
})

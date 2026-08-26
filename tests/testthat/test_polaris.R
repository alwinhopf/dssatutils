test_that("POLARIS tile addressing is stable", {
  expect_equal(dssatutils:::polaris_tile(42.35, -93.40), "lat4243_lon-94-93")
  expect_equal(dssatutils:::polaris_tile(35.91, -101.40), "lat3536_lon-102-101")
  expect_equal(dssatutils:::polaris_tile(40.0, -100.0), "lat4041_lon-100-99")
})

test_that("POLARIS back-transform and water-limit guards are correct", {
  expect_equal(dssatutils:::polaris_backtransform("om", 0), 1)
  expect_equal(dssatutils:::polaris_backtransform("alpha", -2), 0.01)
  expect_equal(dssatutils:::polaris_backtransform("clay", 25), 25)

  wl <- dssatutils:::water_limits(0.08, 0.46, 0.02, 1.4)
  expect_lt(wl$SLLL, wl$SDUL)
  expect_lt(wl$SDUL, wl$SSAT)
  expect_gte(wl$SLLL, 0.02)
  expect_gte(wl$SDUL - wl$SLLL, 0.04 - 1e-10)
  expect_gte(wl$SSAT - wl$SDUL, 0.04 - 1e-10)
})

test_that("POLARIS terra extraction returns raster values, not point IDs", {
  skip_if_not_installed("terra")
  r <- terra::rast(ncols = 2, nrows = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  terra::values(r) <- c(101, 102, 103, 104)
  f <- tempfile(fileext = ".tif")
  terra::writeRaster(r, f, overwrite = TRUE)
  on.exit(unlink(f), add = TRUE)

  pts <- matrix(c(0.5, 1.5,
                  1.5, 0.5), ncol = 2, byrow = TRUE)
  got <- dssatutils:::polaris_extract_values(f, pts, n_expected = 2L, retries = 1L)

  # The old implementation used terra::extract(...)[,1], which can be the ID
  # column (1, 2, ...). Values here deliberately cannot be confused with IDs.
  expect_length(got, 2L)
  expect_true(all(got %in% c(101, 102, 103, 104)))
  expect_false(identical(got, c(1, 2)))
})

test_that("POLARIS failed raster reads terminate and return NA", {
  skip_if_not_installed("terra")
  missing <- tempfile(fileext = ".tif")
  expect_warning(
    got <- dssatutils:::polaris_extract_values(
      missing, matrix(c(0, 0), ncol = 2), n_expected = 1L, retries = 1L
    ),
    "failed after 1 attempt"
  )
  expect_length(got, 1L)
  expect_true(is.na(got))
})

test_that("POLARIS timeout and retry environment settings are validated", {
  old <- Sys.getenv(c("POLARIS_TIMEOUT_SEC", "POLARIS_RETRIES", "POLARIS_PROGRESS_EVERY"), unset = NA)
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, setNames(list(old[[nm]]), nm))
    }
  }, add = TRUE)

  Sys.setenv(POLARIS_TIMEOUT_SEC = "17", POLARIS_RETRIES = "2", POLARIS_PROGRESS_EVERY = "7")
  expect_equal(dssatutils:::polaris_timeout_sec(), 17)
  expect_equal(dssatutils:::polaris_retries(), 2L)
  expect_equal(dssatutils:::polaris_progress_every(), 7L)

  Sys.setenv(POLARIS_TIMEOUT_SEC = "bad", POLARIS_RETRIES = "0", POLARIS_PROGRESS_EVERY = "-1")
  expect_equal(dssatutils:::polaris_timeout_sec(), 45)
  expect_equal(dssatutils:::polaris_retries(), 3L)
  expect_equal(dssatutils:::polaris_progress_every(), 10L)
})

test_that("POLARIS DSSAT writer produces ordered water columns", {
  prof <- data.frame(
    ID = rep("00000001", 2), latitude = rep(42.35, 2), longitude = rep(-93.40, 2),
    depth_bottom = c(5, 15), depth_center = c(2.5, 10),
    SLLL = c(0.177, 0.182), SDUL = c(0.415, 0.421), SSAT = c(0.460, 0.466),
    SSKS = c(1.2, 0.9), bd = c(1.35, 1.38), oc_pct = c(2.1, 1.8),
    clay = c(26, 28), silt = c(41, 40), ph = c(6.2, 6.4)
  )
  d <- tempfile(); dir.create(d); on.exit(unlink(d, recursive = TRUE), add = TRUE)
  dssatutils:::format_dssat_sol_file_polaris(prof, d)
  txt <- readLines(file.path(d, "00000001.SOL"))
  layer <- strsplit(trimws(txt[grepl("^\\s*5\\s", txt)][1]), "\\s+")[[1]]
  expect_lt(as.numeric(layer[3]), as.numeric(layer[4]))
  expect_lt(as.numeric(layer[4]), as.numeric(layer[5]))
})

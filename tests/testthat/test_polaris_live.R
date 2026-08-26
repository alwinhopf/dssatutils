test_that("live POLARIS raster can be sampled through terra", {
  skip_if(Sys.getenv("RUN_POLARIS_LIVE", "false") != "true", "live POLARIS test not enabled")
  skip_if_not_installed("terra")

  old_timeout <- Sys.getenv("POLARIS_TIMEOUT_SEC", unset = NA_character_)
  old_retries <- Sys.getenv("POLARIS_RETRIES", unset = NA_character_)
  on.exit({
    if (is.na(old_timeout)) Sys.unsetenv("POLARIS_TIMEOUT_SEC") else Sys.setenv(POLARIS_TIMEOUT_SEC = old_timeout)
    if (is.na(old_retries)) Sys.unsetenv("POLARIS_RETRIES") else Sys.setenv(POLARIS_RETRIES = old_retries)
  }, add = TRUE)
  Sys.setenv(POLARIS_TIMEOUT_SEC = "30", POLARIS_RETRIES = "2")

  # Hardin County, Iowa: same known-good tile used by the Python POLARIS tests.
  lat <- 42.35
  lon <- -93.40
  tile <- dssatutils:::polaris_tile(lat, lon)
  expect_equal(tile, "lat4243_lon-94-93")

  src <- dssatutils:::polaris_tile_source("clay", "p50", "0_5", tile, cache_dir = NULL)
  expect_match(src, "/vsicurl/")

  value <- suppressWarnings(dssatutils:::polaris_extract_values(
    src,
    matrix(c(lon, lat), ncol = 2),
    n_expected = 1L,
    retries = 2L
  ))

  expect_length(value, 1L)
  expect_true(is.finite(value))
  # Clay is stored linearly as mass percent; broad physical sanity range.
  expect_gte(value, 0)
  expect_lte(value, 100)
})

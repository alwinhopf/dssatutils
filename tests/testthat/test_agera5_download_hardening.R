library(testthat)
library(dssatutils)

test_that("final AgERA5 gridded downloader is self-contained for PSOCK export", {
  body_txt <- paste(deparse(body(dssatutils:::.agera5_download_job)), collapse = "\n")
  expect_match(body_txt, "target = basename\\(job\\$zip_dest\\)")
  expect_false(grepl("basename\\(partial\\)", body_txt))
  expect_false(grepl("\\.agera5_recover_cache_zip\\(", body_txt))
  expect_match(body_txt, "valid_zip <- function")
  expect_match(body_txt, "promote <- function")
})

test_that("AgERA5 time-series downloader requests canonical CSV target", {
  body_txt <- paste(deparse(body(dssatutils:::.agera5_download_timeseries_job)), collapse = "\n")
  expect_match(body_txt, "target = basename\\(dest\\)")
  expect_false(grepl("target = basename\\(partial\\)", body_txt))
  expect_match(body_txt, "paste0\\(dest, \\"\\.csv\\"\\)")
})

test_that("legacy zip.zip cache remains recoverable", {
  skip_if(!nzchar(Sys.which("zip")), "system zip utility is required")
  work <- tempfile("agera5-hardening-")
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  payload <- file.path(work, "payload")
  dir.create(payload)
  writeBin(charToRaw("synthetic netcdf payload"), file.path(payload, "daily.nc"))

  canonical <- file.path(work, "agera5_test.zip")
  legacy <- paste0(canonical, ".zip")
  oldwd <- setwd(payload)
  on.exit(setwd(oldwd), add = TRUE)
  suppressWarnings(utils::zip(zipfile = legacy, files = "daily.nc", flags = "-j"))
  setwd(oldwd)

  expect_true(dssatutils:::.agera5_valid_zip(legacy))
  recovered <- dssatutils:::.agera5_recover_cache_zip(canonical)
  expect_equal(recovered, canonical)
  expect_true(dssatutils:::.agera5_valid_zip(canonical))
  expect_false(file.exists(legacy))
})

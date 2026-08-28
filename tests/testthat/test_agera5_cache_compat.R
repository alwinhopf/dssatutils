library(testthat)
library(dssatutils)

test_that("AgERA5 recovers legacy ecmwfr .zip.zip cache artifacts", {
  skip_if(!nzchar(Sys.which("zip")), "system zip utility is required for this archive test")

  work <- tempfile("agera5-cache-")
  dir.create(work, recursive = TRUE)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)

  payload_dir <- file.path(work, "payload")
  dir.create(payload_dir)
  writeBin(charToRaw("synthetic netcdf payload"), file.path(payload_dir, "daily.nc"))

  canonical <- file.path(work, "agera5_test.zip")
  legacy <- paste0(canonical, ".zip")
  oldwd <- setwd(payload_dir)
  on.exit(setwd(oldwd), add = TRUE)
  suppressWarnings(utils::zip(zipfile = legacy, files = "daily.nc", flags = "-j"))
  setwd(oldwd)

  expect_true(file.exists(legacy))
  expect_false(file.exists(canonical))
  expect_true(dssatutils:::.agera5_valid_zip(legacy))

  recovered <- dssatutils:::.agera5_recover_cache_zip(canonical)
  expect_equal(recovered, canonical)
  expect_true(file.exists(canonical))
  expect_false(file.exists(legacy))
  expect_true(dssatutils:::.agera5_valid_zip(canonical))
})

test_that("AgERA5 downloader requests the canonical zip target", {
  downloader_body <- paste(deparse(body(dssatutils:::.agera5_download_job)), collapse = "\n")
  expect_match(downloader_body, "target = basename\\(job\\$zip_dest\\)")
  expect_false(grepl("basename\\(partial\\)", downloader_body))
})

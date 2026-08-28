library(testthat)
library(dssatutils)

test_that("final AgERA5 gridded downloader is self-contained for PSOCK export", {
  body_txt <- paste(deparse(body(dssatutils:::.agera5_download_job)), collapse = "\n")
  expect_match(body_txt, "\\.agera5_gridded_request\\(job, basename\\(job\\$zip_dest\\)\\)")
  expect_false(grepl("basename\\(partial\\)", body_txt))
  expect_false(grepl("\\.agera5_recover_cache_zip\\(", body_txt))
  expect_match(body_txt, "valid_zip <- function")
  expect_match(body_txt, "promote <- function")
})

test_that("AgERA5 time-series downloader requests canonical CSV target", {
  body_txt <- paste(deparse(body(dssatutils:::.agera5_download_timeseries_job)), collapse = "\n")
  expect_match(body_txt, "\\.agera5_timeseries_request\\(job, basename\\(dest\\)\\)")
  expect_false(grepl("target = basename\\(partial\\)", body_txt))
  expect_match(body_txt, "paste0\\(dest")
  expect_match(body_txt, "\\.csv")
})

test_that("AgERA5 request builders match current CDS catalogue forms", {
  gridded_job <- dssatutils:::.agera5_job(
    "TMAX", 1993, dssatutils:::.agera5_vars$TMAX,
    c(33.8816, -102.7220, 33.4816, -102.3220), tempdir()
  )
  gridded <- dssatutils:::.agera5_gridded_request(gridded_job, "agera5.zip")
  expect_equal(gridded$dataset_short_name, "sis-agrometeorological-indicators")
  expect_equal(gridded$version, "2_0")
  expect_equal(gridded$variable, "2m_temperature")
  expect_equal(gridded$statistic, "24_hour_maximum")
  expect_equal(length(gridded$month), 12)
  expect_equal(length(gridded$day), 31)

  timeseries_job <- list(
    year = 1993, area = c(33.8816, -102.7220, 33.4816, -102.3220),
    data_format = "csv"
  )
  timeseries <- dssatutils:::.agera5_timeseries_request(timeseries_job, "agera5.csv")
  expect_equal(timeseries$dataset_short_name,
               "sis-agrometeorological-indicators-timeseries")
  expect_equal(timeseries$data_format, "csv")
  expect_equal(timeseries$date, c("1993-01-01", "1993-12-31"))
  expect_equal(length(timeseries$variable), 7)
  expect_false("version" %in% names(timeseries))
})
test_that("AgERA5 retries a throttled submission with server-aware backoff", {
  calls <- 0L
  waits <- numeric()
  request_fn <- function(request, path) {
    calls <<- calls + 1L
    if (calls == 1L) stop("429 Rate limit exceeded. Please wait 1 seconds.")
    "download.zip"
  }

  result <- dssatutils:::.agera5_wf_request_with_retry(
    request = list(), path = tempdir(), target = "download.zip",
    request_fn = request_fn,
    resume_fn = function(...) stop("resume should not be called"),
    sleep_fn = function(seconds) waits <<- c(waits, seconds),
    max_attempts = 3L
  )

  expect_equal(result$value, "download.zip")
  expect_null(result$error)
  expect_equal(calls, 2L)
  expect_equal(waits, 2)
})

test_that("AgERA5 resumes the submitted CDS job instead of duplicating it after 429", {
  job_url <- paste0(
    "https://cds.climate.copernicus.eu/api/retrieve/v1/jobs/",
    "aa31ed9e-6a7f-423d-a9e6-a82b35d3f252"
  )
  submissions <- 0L
  resumes <- 0L
  request_fn <- function(request, path) {
    submissions <<- submissions + 1L
    stop(paste("429 Rate limit exceeded. Please wait 0 seconds.", job_url))
  }
  resume_fn <- function(url, path, target) {
    resumes <<- resumes + 1L
    expect_equal(url, job_url)
    if (resumes == 1L) stop("Your requested file is unavailable - check url")
    invisible(list(asset = target))
  }

  result <- dssatutils:::.agera5_wf_request_with_retry(
    request = list(), path = tempdir(), target = "download.zip",
    request_fn = request_fn, resume_fn = resume_fn,
    sleep_fn = function(seconds) NULL, max_attempts = 4L
  )

  expect_null(result$error)
  expect_equal(result$job_url, job_url)
  expect_equal(submissions, 1L)
  expect_equal(resumes, 2L)
})

test_that("AgERA5 does not retry permanent request errors", {
  calls <- 0L
  result <- dssatutils:::.agera5_wf_request_with_retry(
    request = list(), path = tempdir(), target = "download.zip",
    request_fn = function(request, path) {
      calls <<- calls + 1L
      stop("required licences not accepted")
    },
    sleep_fn = function(seconds) stop("sleep should not be called"),
    max_attempts = 3L
  )
  expect_equal(calls, 1L)
  expect_match(result$error, "licences not accepted")
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

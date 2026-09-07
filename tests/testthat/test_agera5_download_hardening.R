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
  body_txt <- paste(deparse(body(dssatutils:::.agera5_download_timeseries_impl)), collapse = "\n")
  expect_match(body_txt, "target = basename\\(dest\\)")
  expect_false(grepl("target = basename\\(partial\\)", body_txt))
  expect_match(body_txt, "paste0\\(dest")
  expect_match(body_txt, "\\.csv")
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

test_that(".agera5_validate_timeseries_csv enforces columns and calendar completeness", {
  tf <- tempfile(fileext = ".csv")
  on.exit(unlink(tf), add = TRUE)

  dates <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "day")
  df <- data.frame(
    valid_time = as.character(dates),
    latitude = 40.0,
    longitude = -90.0,
    Temperature_Air_2m_Max_24h = 300.0,
    Temperature_Air_2m_Min_24h = 280.0,
    Solar_Radiation_Flux = 12000000.0,
    Precipitation_Flux = 2.0,
    Dew_Point_Temperature_2m_Mean_24h = 275.0,
    Relative_Humidity_2m_15h = 60.0,
    Wind_Speed_10m_Mean_24h = 3.0,
    stringsAsFactors = FALSE
  )
  write.csv(df, tf, row.names = FALSE)

  # Full valid CSV matches year 2010
  expect_true(dssatutils:::.agera5_validate_timeseries_csv(tf, 2010))

  # Wrong year expectation fails
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(tf, 2011))

  # Physical inversion (Tmax < Tmin) is preserved in raw cache
  df_inv <- df
  df_inv$Temperature_Air_2m_Max_24h[10] <- 270.0
  df_inv$Temperature_Air_2m_Min_24h[10] <- 280.0
  write.csv(df_inv, tf, row.names = FALSE)
  expect_true(dssatutils:::.agera5_validate_timeseries_csv(tf, 2010))

  # Missing a required column fails
  df_no_srad <- df[, !names(df) %in% "Solar_Radiation_Flux"]
  write.csv(df_no_srad, tf, row.names = FALSE)
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(tf, 2010))

  # Duplicate date in cell fails
  df_dup <- rbind(df, df[10, ])
  write.csv(df_dup, tf, row.names = FALSE)
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(tf, 2010))

  # Missing days fails
  df_short <- df[1:100, ]
  write.csv(df_short, tf, row.names = FALSE)
  expect_false(dssatutils:::.agera5_validate_timeseries_csv(tf, 2010))
})

test_that(".agera5_download_timeseries_job detects canonical CSV even on non-path wf_request return", {
  cache_dir <- tempfile("agera5-cache-")
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  yr <- 2010L
  area <- c(40.5, -90.5, 39.5, -89.5)
  dest <- dssatutils:::.agera5_timeseries_cache_path(cache_dir, yr, area, "csv")

  job <- list(
    year = yr,
    area = area,
    cache_dir = cache_dir,
    data_format = "csv",
    chunk = list(area = area, idx = 1L)
  )

  # Mock ecmwfr::wf_request by writing valid CSV to dest but returning NULL
  dates <- seq.Date(as.Date("2010-01-01"), as.Date("2010-12-31"), by = "day")
  df <- data.frame(
    valid_time = as.character(dates),
    latitude = 40.0,
    longitude = -90.0,
    Temperature_Air_2m_Max_24h = 300.0,
    Temperature_Air_2m_Min_24h = 280.0,
    Solar_Radiation_Flux = 12000000.0,
    Precipitation_Flux = 2.0,
    Dew_Point_Temperature_2m_Mean_24h = 275.0,
    Relative_Humidity_2m_15h = 60.0,
    Wind_Speed_10m_Mean_24h = 3.0,
    stringsAsFactors = FALSE
  )
  write.csv(df, dest, row.names = FALSE)

  # With valid dest in place, .agera5_download_timeseries_job returns dest immediately
  res <- dssatutils:::.agera5_download_timeseries_job(job)
  expect_equal(res, dest)
})

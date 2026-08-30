library(testthat)
library(dssatutils)

write_sample_wth <- function(path, rows) {
  writeLines(c(
    "$WEATHER DATA: test",
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    "  TEST   0.0000   0.0000   -99  10.0  20.0   2.0  10.0",
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    rows
  ), path)
}

test_that("fixed-width adjacent negative weather values are valid", {
  path <- tempfile(fileext = ".WTH")
  rows <- c(
    sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f", "2024001", 12, -10.2, -12.3, 0, -15, 40, 3),
    sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f", "2024002", 13, -9, -11, 0, -14, 42, 3.2)
  )
  write_sample_wth(path, rows)
  expect_true(is_wth_valid(path, end_year = 2024))
})

test_that("weather validation rejects date gaps", {
  path <- tempfile(fileext = ".WTH")
  write_sample_wth(path, c(
    "2024001 12.0 10.0 1.0 0.0 0.0 40.0 3.0",
    "2024003 12.0 10.0 1.0 0.0 0.0 40.0 3.0"
  ))
  expect_false(is_wth_valid(path, end_year = 2024))
})

test_that("weather validation rejects absolute-zero temperatures", {
  path <- tempfile(fileext = ".WTH")
  row <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                 "2018001", 0, -273.1, -273.1, 0, -273.1, 0, 0)
  write_sample_wth(path, row)
  expect_false(is_wth_valid(path, end_year = 2018))
})

test_that("weather validation can require complete core forcing", {
  path <- tempfile(fileext = ".WTH")
  write_sample_wth(path, c(
    "2024001 -99 10.0 1.0 0.0 0.0 40.0 -99",
    "2024002 12.0 10.0 1.0 0.0 0.0 40.0 -99"
  ))
  expect_true(is_wth_valid(path, end_year = 2024))
  expect_false(is_wth_valid(
    path,
    end_year = 2024,
    required_columns = c("SRAD", "TMAX", "TMIN", "RAIN")
  ))
})

test_that("weather validation can require all AgERA5 forcing", {
  path <- tempfile(fileext = ".WTH")
  write_sample_wth(path, c(
    "2024001 12.0 10.0 1.0 0.0 -99 -99 -99",
    "2024002 12.0 10.0 1.0 0.0 -99 -99 -99"
  ))
  core <- c("SRAD", "TMAX", "TMIN", "RAIN")
  agera5 <- c(core, "TDEW", "RH2M", "WIND")

  expect_true(is_wth_valid(path, end_year = 2024, required_columns = core))
  expect_false(is_wth_valid(path, end_year = 2024, required_columns = agera5))
})

test_that("AgERA5 writer defers physical validation to the shared validator", {
  wd <- data.frame(
    DATE = c("2018001", "2018002"), YEAR = c(2018, 2018), MM = c(1, 1),
    SRAD = 12, TMAX = c(5.8, 10), TMIN = c(6, 2), RAIN = 0,
    TDEW = 0, RH2M = 60, WIND = 3
  )
  writer <- getFromNamespace(".agera5_write_wth", "dssatutils")
  path <- tempfile(fileext = ".WTH")
  out_dir <- dirname(path)
  generated <- writer(wd, tools::file_path_sans_ext(basename(path)), 33.7, -102.5, out_dir)

  expect_true(file.exists(generated))
  expect_true(any(grepl("   5.8   6.0", readLines(generated), fixed = TRUE)))
  expect_false(is_wth_valid(
    generated,
    end_year = 2018,
    required_columns = c("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
  ))
})

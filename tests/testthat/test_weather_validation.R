library(testthat)


write_sample_wth <- function(path, rows) {
  writeLines(c(
    "$WEATHER DATA: test",
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    "  TEST   0.0000   0.0000   -99  10.0  20.0   2.0  10.0",
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    rows
  ), path)
}

test_that("wind validation uses DSSAT kilometres per day", {
  path <- tempfile(fileext = ".WTH")
  write_sample_wth(path, "2024001 12.0 10.0 1.0 0.0 0.0 40.0 1000.0")
  expect_true(is_wth_valid(path, end_date = "2024-01-01"))
  write_sample_wth(path, "2024001 12.0 10.0 1.0 0.0 0.0 40.0 9000.0")
  expect_false(is_wth_valid(path, end_date = "2024-01-01"))
})

test_that("fixed-width adjacent negative weather values are valid", {
  path <- tempfile(fileext = ".WTH")
  rows <- c(
    sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f", "2024001", 12, -10.2, -12.3, 0, -15, 40, 3),
    sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f", "2024002", 13, -9, -11, 0, -14, 42, 3.2)
  )
  write_sample_wth(path, rows)
  expect_true(is_wth_valid(path, end_date = "2024-01-02"))
})

test_that("weather validation rejects date gaps", {
  path <- tempfile(fileext = ".WTH")
  write_sample_wth(path, c(
    "2024001 12.0 10.0 1.0 0.0 0.0 40.0 3.0",
    "2024003 12.0 10.0 1.0 0.0 0.0 40.0 3.0"
  ))
  expect_false(is_wth_valid(path, end_date = "2024-01-02"))
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
  expect_true(is_wth_valid(path, end_date = "2024-01-02"))
  expect_false(is_wth_valid(
    path,
    end_date = "2024-01-02",
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

  expect_true(is_wth_valid(path, end_date = "2024-01-02", required_columns = core))
  expect_false(is_wth_valid(path, end_date = "2024-01-02", required_columns = agera5))
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


test_that("AgERA5 writes wind run and preserves missing sentinel", {
  wd <- data.frame(DATE = c("2018001", "2018002"), YEAR = 2018, MM = 1,
                   SRAD = 12, TMAX = 10, TMIN = 2, RAIN = 0,
                   TDEW = 0, RH2M = 60, WIND = c(3, -99))
  writer <- getFromNamespace(".agera5_write_wth", "dssatutils")
  generated <- writer(wd, "WIND", 30, -87, tempdir())
  lines <- readLines(generated)
  lines <- lines[grepl("^[0-9]{7}", lines)]
  expect_equal(as.numeric(substr(lines, 44, 49)), c(259.2, -99))
  expect_true(is_wth_valid(generated, end_date = "2018-01-02"))
})

test_that("AgERA5 preserves nonzero temperatures near freezing", {
  wd <- data.frame(DATE = c("2007361", "2007362", "2007363"), YEAR = 2007, MM = 12,
                   SRAD = 3.3, TMAX = c(0.0432, 0, 2.12), TMIN = c(0.03686, 0, -1.21),
                   RAIN = 2.2, TDEW = -3.7, RH2M = 76.4, WIND = 2.8)
  writer <- getFromNamespace(".agera5_write_wth", "dssatutils")
  lines <- tail(readLines(writer(wd, "PREC", 36.52579, -98.254504, tempdir())), 3)
  expect_equal(nchar(lines), rep(49L, 3))
  expect_equal(as.numeric(substr(lines, 14, 19)), c(0.04, 0, 2.1))
  expect_equal(as.numeric(substr(lines, 20, 25)), c(0.04, 0, -1.2))
})

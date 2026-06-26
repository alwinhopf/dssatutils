test_that("temperature inversion repair uses neighboring Tmax/Tmin means", {
  wth <- tempfile(fileext = ".WTH")
  log_file <- tempfile(fileext = ".log")
  writeLines(c(
    "$WEATHER DATA: TEST",
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    "  TEST   0.0000   0.0000   -99  20.0  10.0   2.0   2.0",
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    "2024001  15.0  20.0  10.0   0.0   8.0  60.0   2.0",
    "2024002  15.0  22.0  12.0   0.0   9.0  60.0   2.0",
    "2024003  15.0   5.0  15.0   0.0  10.0  60.0   2.0",
    "2024004  15.0  24.0  14.0   0.0  11.0  60.0   2.0",
    "2024005  15.0  26.0  16.0   0.0  12.0  60.0   2.0",
    "2024006  15.0  28.0  18.0   0.0  13.0  60.0   2.0"
  ), wth, useBytes = TRUE)

  summary <- repair_weather_file_temperature_inversions(
    wth,
    max_gap_days = 3,
    window_days = 2,
    log_file = log_file
  )

  expect_equal(summary$status, "repaired")
  expect_equal(summary$repaired_count, 1L)

  rows <- readLines(wth, warn = FALSE)
  rows <- rows[grepl("^\\s*2024", rows)]
  dat <- utils::read.table(
    text = paste(rows, collapse = "\n"),
    col.names = c("DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"),
    colClasses = c("character", rep("numeric", 7)),
    stringsAsFactors = FALSE
  )
  repaired <- dat[dat$DATE == "2024003", ]

  expect_equal(repaired$TMAX, 23)
  expect_equal(repaired$TMIN, 13)
  expect_lte(repaired$TMIN, repaired$TMAX)
  expect_true(any(grepl("issue=TMIN_GT_TMAX status=repaired", readLines(log_file))))
})

test_that("date gap repair inserts missing DSSAT date row from neighbor means", {
  wth <- tempfile(fileext = ".WTH")
  log_file <- tempfile(fileext = ".log")
  rows <- c(
    "$WEATHER DATA: TEST",
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    "  TEST   0.0000   0.0000   -99  20.0  10.0   2.0   2.0",
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    "2024001  15.0  20.0  10.0   0.0   8.0  60.0   2.0",
    "2024002  15.0  22.0  12.0   0.0   9.0  60.0   2.0",
    "2024004  15.0  24.0  14.0   0.0  11.0  60.0   2.0",
    "2024005  15.0  26.0  16.0   0.0  12.0  60.0   2.0",
    "2024006  15.0  28.0  18.0   0.0  13.0  60.0   2.0"
  )
  writeLines(rows, wth, useBytes = TRUE)

  summary <- repair_weather_file_date_gaps(
    wth,
    max_gap_days = 3,
    window_days = 2,
    log_file = log_file
  )

  expect_equal(summary$status, "repaired")
  expect_equal(summary$repaired_count, 1L)

  rows <- readLines(wth, warn = FALSE)
  rows <- rows[grepl("^\\s*2024", rows)]
  dat <- utils::read.table(
    text = paste(rows, collapse = "\n"),
    col.names = c("DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"),
    colClasses = c("character", rep("numeric", 7)),
    stringsAsFactors = FALSE
  )
  inserted <- dat[dat$DATE == "2024003", ]

  expect_equal(inserted$TMAX, 23)
  expect_equal(inserted$TMIN, 13)
  expect_true(any(grepl("issue=DATE_GAP status=repaired", readLines(log_file))))
})

test_that("weather quality audit flags suspicious rows without modifying file", {
  wth <- tempfile(fileext = ".WTH")
  writeLines(c(
    "$WEATHER DATA: TEST",
    "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
    "  TEST   0.0000   0.0000   -99  20.0  10.0   2.0   2.0",
    "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
    "2024001  15.0  20.0  10.0   0.0   8.0  60.0   2.0",
    "2024002  15.0  22.0  12.0   0.0   9.0  60.0   2.0",
    "2024003  15.0   5.0  15.0   0.0  10.0  60.0   2.0",
    "2024004  15.0  24.0  14.0   0.0  11.0  60.0   2.0"
  ), wth, useBytes = TRUE)

  audit <- audit_weather_file_quality(wth, flatline_days = 3)

  expect_true("tmin_gt_tmax" %in% audit$issue)
  expect_true("RAIN_flatline" %in% audit$issue)
})

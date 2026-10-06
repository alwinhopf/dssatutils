library(testthat)
library(dssatutils)

test_that("FAO-56 Example 8 extraterrestrial radiation verification", {
  res <- extraterrestrial_radiation(-20, 246)
  expect_equal(res$Ra, 32.2, tolerance = 0.1)
  expect_equal(round(res$Ra, 2), 32.19)
  expect_true(res$ws > 0 && res$ws < pi)
})

test_that("Extraterrestrial radiation handles polar limits safely", {
  # North Pole winter solstice (24h darkness, Ra = 0)
  np_win <- extraterrestrial_radiation(90, 355)
  expect_equal(np_win$Ra, 0)
  expect_equal(np_win$ws, 0)

  # North Pole summer solstice (24h sunlight, ws = pi, Ra > 0)
  np_sum <- extraterrestrial_radiation(90, 172)
  expect_gt(np_sum$Ra, 0)
  expect_equal(np_sum$ws, pi, tolerance = 1e-6)

  # South Pole winter solstice (24h darkness, Ra = 0)
  sp_win <- extraterrestrial_radiation(-90, 172)
  expect_equal(sp_win$Ra, 0)
  expect_equal(sp_win$ws, 0)

  # South Pole summer solstice (24h sunlight, ws = pi, Ra > 0)
  sp_sum <- extraterrestrial_radiation(-90, 355)
  expect_gt(sp_sum$Ra, 0)
  expect_equal(sp_sum$ws, pi, tolerance = 1e-6)
})

test_that("Bristow-Campbell satisfies physical bounds 0 <= Rs <= 0.80 * Ra", {
  dates <- seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  doy <- as.integer(format(dates, "%j"))
  tmax <- 20 + 10 * sin(2 * pi * doy / 365)
  tmin <- tmax - seq(4, 20, length.out = length(dates))

  for (lat in c(-45, -15, 0, 30, 50)) {
    rs <- estimate_srad_bristow_campbell(dates, tmax, tmin, lat)
    e <- extraterrestrial_radiation(lat, doy)
    expect_true(all(is.finite(rs)))
    expect_true(all(rs >= 0))
    expect_true(all(rs <= 0.80 * e$Ra + 1e-6))
  }
})

test_that("Bristow-Campbell responds monotonically to diurnal range", {
  dates <- seq(as.Date("2020-07-01"), by = "day", length.out = 10)
  lat <- 35.0

  rs_low <- estimate_srad_bristow_campbell(dates, rep(22, 10), rep(20, 10), lat)
  rs_high <- estimate_srad_bristow_campbell(dates, rep(35, 10), rep(17, 10), lat)

  expect_true(all(rs_high > rs_low))
})

test_that("Bristow-Campbell handles date gaps and missing temperatures", {
  dates <- as.Date(c("2020-01-01", "2020-01-02", "2020-01-10", "2020-01-11"))
  tmax <- c(25, NA, 22, 24)
  tmin <- c(10, 12, 10, 11)

  rs <- estimate_srad_bristow_campbell(dates, tmax, tmin, 30.0)
  expect_equal(length(rs), 4)
  expect_true(is.finite(rs[1]) && rs[1] > 0)
  expect_true(is.na(rs[2]))
  expect_true(is.finite(rs[3]) && rs[3] > 0)
  expect_true(is.finite(rs[4]) && rs[4] > 0)
})

test_that("Bristow-Campbell cross-language numerical parity with Python", {
  # Deterministic test sequence
  dates <- as.Date("2020-01-01") + 0:4
  tmax <- c(15.0, 18.0, NA, 20.0, 22.0)
  tmin <- c(5.0, 7.0, 6.0, NA, 10.0)
  rs <- estimate_srad_bristow_campbell(dates, tmax, tmin, 30.0)

  # Python calculated: [11.00180729, 13.64421333, nan, nan, 14.08366375]
  expect_equal(rs[1], 11.00180729, tolerance = 1e-5)
  expect_equal(rs[2], 13.64421333, tolerance = 1e-5)
  expect_true(is.na(rs[3]))
  expect_true(is.na(rs[4]))
  expect_equal(rs[5], 14.08366375, tolerance = 1e-5)
})

test_that("PRISM srad_method parameter and WTH headers", {
  pts <- sf::st_sf(
    ID = "P1",
    geometry = sf::st_sfc(sf::st_point(c(-85.0, 32.5)), crs = 4326)
  )

  # Invalid method throws error
  expect_error(
    process_weather_prism(pts, 2020, 2020, tempfile(), "ID", "LAT", "LONG", 1,
                          tempfile(), tempfile(), srad_method = "invalid_method"),
    "Unknown srad_method"
  )
})

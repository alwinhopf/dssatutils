library(testthat)
library(dssatutils)

# Offline tests for the new gNATSGO (soil), DWD and E-OBS (weather) sources.
# The network-bound lookups are isolated from the physics/formatting helpers,
# which are exercised here with deterministic synthetic inputs (no network, no
# DSSAT install). Mirrors python/tests/test_new_sources.py.

# ===================================================================
# DWD: solar radiation from sunshine duration + dew point
# ===================================================================
test_that("DWD extraterrestrial radiation peaks in summer", {
  Ra_win <- dssatutils:::.dwd_extraterrestrial(50, 15)$Ra
  Ra_sum <- dssatutils:::.dwd_extraterrestrial(50, 196)$Ra
  expect_gt(Ra_sum, Ra_win)
  expect_gt(Ra_win, 0)
})

test_that("DWD SRAD from sunshine is in a physical range", {
  e <- dssatutils:::.dwd_extraterrestrial(50, 196)
  N <- 24 / pi * e$ws
  rs_full <- dssatutils:::.dwd_srad_from_sunshine(50, 196, N)       # clear sky
  expect_true(rs_full > 18 && rs_full < 32)
  rs_zero <- dssatutils:::.dwd_srad_from_sunshine(50, 196, 0)       # diffuse only
  expect_true(rs_zero > 0 && rs_zero < rs_full)
  expect_true(is.na(dssatutils:::.dwd_srad_from_sunshine(50, 196, NA)))
})

test_that("DWD dew point from vapour pressure", {
  td <- dssatutils:::.dwd_tdew_from_vp(12.27)
  expect_lt(abs(td - 10.0), 0.5)
})

# ===================================================================
# E-OBS: dew point from RH + .WTH writer
# ===================================================================
test_that("E-OBS dew point from RH and writer", {
  td <- dssatutils:::.eobs_tdew_from_rh(23.0, 65.0)
  expect_true(td > 15 && td < 17)

  dates <- seq(as.Date("2019-01-01"), as.Date("2019-12-31"), by = "day")
  df <- data.frame(
    DATE = sprintf("%d%03d", as.integer(format(dates, "%Y")), as.integer(format(dates, "%j"))),
    YEAR = as.integer(format(dates, "%Y")), MM = as.integer(format(dates, "%m")),
    SRAD = 15, TMAX = 20, TMIN = 9, RAIN = 1, TDEW = 8, RH2M = 70, WIND = 3)
  work <- tempfile(); dir.create(work); on.exit(unlink(work, recursive = TRUE))
  dssatutils:::.eobs_write_wth(df, "EOBS_T", 50.0, 8.0, work)
  out <- file.path(work, "EOBS_T.WTH")
  expect_true(file.exists(out))
  ln <- readLines(out)
  expect_true(grepl("\\$WEATHER DATA: E-OBS", ln[1]))
  expect_true(any(grepl("@  DATE", ln)))
  expect_equal(sum(grepl("^2019", ln)), 365)
})

# ===================================================================
# gNATSGO: the .SOL writer (no network)
# ===================================================================
test_that("gNATSGO .SOL writer produces a valid profile", {
  `%>%` <- dplyr::`%>%`
  profile <- data.frame(
    ID = "00000099", latitude = 42.0, longitude = -93.0,
    depth_range = c("0-5cm", "5-20cm"),
    clay_pct = c(25, 26), sand_pct = c(35, 34), silt_pct = c(40, 40),
    om_pct = c(3.0, 2.5), bulk_density = c(1.3, 1.35),
    SLLL = c(0.12, 0.13), SDUL = c(0.28, 0.29), SSAT = c(0.45, 0.45))
  work <- tempfile(); dir.create(work); on.exit(unlink(work, recursive = TRUE))
  dssatutils:::format_dssat_soil_gnatsgo(profile, work)
  out <- file.path(work, "00000099.SOL")
  expect_true(file.exists(out))
  txt <- readLines(out)
  expect_true(any(grepl("USA gNATSGO Soil Profiles", txt)))
  expect_true(any(grepl("@  SLB", txt)))
})

# ===================================================================
# iSDAsoil (Africa): uint8 back-transformation
# ===================================================================
test_that("iSDAsoil back-transform converts stored values correctly", {
  expect_equal(dssatutils:::.isda_back_transform("clay_content", 42), 42)
  expect_equal(dssatutils:::.isda_back_transform("bulk_density", 117), 1.17)
  expect_equal(dssatutils:::.isda_back_transform("carbon_organic", 30), exp(3) - 1, tolerance = 1e-6)
  expect_true(is.na(dssatutils:::.isda_back_transform("clay_content", 255)))  # nodata
})

# ===================================================================
# Xavier (Brazil): Rs passes through (already MJ/m^2/day)
# ===================================================================
test_that("Xavier writer does not rescale solar radiation", {
  dates <- seq(as.Date("2015-01-01"), as.Date("2015-12-31"), by = "day")
  df <- data.frame(
    DATE = sprintf("%d%03d", as.integer(format(dates, "%Y")), as.integer(format(dates, "%j"))),
    YEAR = as.integer(format(dates, "%Y")), MM = as.integer(format(dates, "%m")),
    SRAD = 22, TMAX = 32, TMIN = 22, RAIN = 3, TDEW = 20, RH2M = 75, WIND = 2.5)
  work <- tempfile(); dir.create(work); on.exit(unlink(work, recursive = TRUE))
  dssatutils:::.xavier_write_wth(df, "BR1", -15.8, -47.9, work)
  ln <- readLines(file.path(work, "BR1.WTH"))
  expect_true(grepl("BR-DWGD/Xavier", ln[1]))
  expect_equal(as.numeric(substr(ln[5], 8, 13)), 22.0, tolerance = 0.05)
})

# ===================================================================
# CMFD (China): RH/TDEW from specific humidity
# ===================================================================
test_that("CMFD derives plausible RH/TDEW from specific humidity", {
  rt <- dssatutils:::.cmfd_rh_tdew(0.004, 10.0, 90000.0)
  expect_true(rt$rh > 30 && rt$rh < 70)
  expect_true(rt$tdew < 10)
})

# ===================================================================
# LUCAS (Europe): .SOL writer
# ===================================================================
test_that("LUCAS .SOL writer flags topsoil extrapolation", {
  profile <- data.frame(
    ID = "L1", latitude = 50.8, longitude = 6.1, depth_top = c(0, 20), depth_bottom = c(20, 150),
    clay_pct = 22, sand_pct = 40, silt_pct = 38, om_pct = 3.1, bulk_density = 1.5,
    SLLL = 0.13, SDUL = 0.26, SSAT = 0.40)
  work <- tempfile(); dir.create(work); on.exit(unlink(work, recursive = TRUE))
  dssatutils:::format_dssat_soil_lucas(profile, work)
  txt <- readLines(file.path(work, "L1.SOL"))
  expect_true(any(grepl("LUCAS Topsoil", txt)))
  expect_true(any(grepl("EXTRAPOLATED", txt)))
})

# ===================================================================
# Public exports present
# ===================================================================
test_that("new public entry points are exported", {
  for (fn in c("process_soils_gnatsgo", "process_weather_dwd", "process_weather_eobs",
               "process_soils_isdasoil", "process_soils_lucas",
               "process_weather_xavier", "process_weather_cmfd")) {
    expect_true(exists(fn, where = asNamespace("dssatutils")))
  }
})

test_that("declared units give consistent daily energy and pressure", {
  cv <- dssatutils:::weather_convert_units
  expect_equal(cv(20, "MJ m-2 day-1", "srad"), 20)
  expect_equal(cv(20000, "kJ m-2 day-1", "srad"), 20)
  expect_equal(cv(20000000, "J m-2", "srad"), 20)
  expect_equal(cv(1200, "Pa", "vp"), cv(12, "hPa", "vp"))
  expect_equal(cv(1.2, "kPa", "vp"), cv(12, "hPa", "vp"))
  expect_equal(cv(3.6, "km/h", "wind", wind_height_m = 2), 1)
  expect_error(cv(1, "m/s", "wind"), "height_m")
  expect_error(cv(20, "unknown", "srad"), "units")
})

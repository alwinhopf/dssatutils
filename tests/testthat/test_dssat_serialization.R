# These fixtures and invariants mirror test_dssat_serialization.py.
native_values <- function(header, line) {
  positions <- gregexpr("[^[:space:]]+", header)[[1]]
  lengths <- attr(positions, "match.length")
  ends <- positions + lengths - 1L
  names <- regmatches(header, gregexpr("[^[:space:]]+", header))[[1]]
  idx <- seq.int(2L, length(names))
  values <- suppressWarnings(as.numeric(substring(line, ends[idx - 1L] + 2L, ends[idx])))
  stats::setNames(values, names[idx])
}

test_that("numeric fields preserve the native separator and rounding", {
  fixture <- read.csv(testthat::test_path("..", "fixtures", "dssat_numeric_fields.csv"),
                      colClasses = "character", na.strings = NULL)
  fields <- .format_wth_value(fixture$value, as.integer(fixture$decimals))
  expect_true(all(nchar(fields) == 6L & substr(fields, 1, 1) == " "))
  expect_equal(trimws(fields), fixture$expected)
  expect_error(.format_wth_value(100000), "cannot fit")
  expect_equal(.weather_wind_run(c(3, -99, NA, Inf)), c(259.2, -99, -99, -99))
})

test_that("provider wind units and daily columns match native DSSAT", {
  df <- data.frame(DATE=c("2024001", "2024002"), YEAR=2024, MM=1, DOY=1:2,
                   SRAD=12, TMAX=-5, TMIN=-12.3, RAIN=c(999.96, 1200),
                   TDEW=-13, RH2M=76.9, WIND=c(13, -99))
  out <- tempfile(); dir.create(out)
  on.exit(unlink(out, recursive=TRUE))
  writers <- list(.cmfd_write_wth, .eobs_write_wth, .xavier_write_wth, .agera5_write_wth)
  writers <- c(writers, list(
    function(df, pid, lat, lon, out) .dwd_write_wth(df, pid, lat, lon, 100, out),
    function(df, pid, lat, lon, out) weather_write_wth(df, pid, lat, lon, out, "test", "TEST"),
    function(df, pid, lat, lon, out) .write_dssat_weather_file(df, lat, lon, file.path(out, "TEST.WTH"), pid)))
  for (writer in writers) {
    writer(df, "TEST", 50, 10, out)
    lines <- readLines(file.path(out, "TEST.WTH"))
    header <- lines[grepl("^@  DATE", lines)]
    rows <- lines[grepl("^[0-9]{7}", lines)]
    expect_equal(nchar(rows), c(49L, 49L))
    values <- lapply(rows, function(line) native_values(header, line))
    expect_equal(vapply(values, function(v) v[["WIND"]], numeric(1)), c(1123, -99))
    expect_equal(vapply(values, function(v) v[["RAIN"]], numeric(1)), c(1000, 1200))
    expect_equal(vapply(values, function(v) v[["TMIN"]], numeric(1)), c(-12.3, -12.3))
  }
})

test_that("soil conductivity fits the native field at the rounding boundary", {
  profile <- data.frame(ID="TEST", latitude=50, longitude=10, depth_bottom=c(20,50), depth_range=c("0-20cm", "20-50cm"),
                        SLLL=.1, SDUL=.2, SSAT=.4, SSKS=c("99.996","150"),
                        bulk_density=1.4, om_pct=2, clay_pct=20, silt_pct=30)
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive=TRUE))
  for (writer in list(format_dssat_soil_isdasoil, format_dssat_soil_lucas, format_dssat_soil_single, format_dssat_soil_gnatsgo)) {
    path <- file.path(out, "TEST.SOL"); unlink(path)
    writer(profile, out)
    lines <- readLines(path); index <- grep("^@  SLB", lines)
    rows <- lines[seq.int(index + 1L, length(lines))]; rows <- rows[nzchar(trimws(rows))]
    values <- lapply(rows, function(line) native_values(lines[index], line))
    expect_equal(vapply(values, function(v) v[["SSKS"]], numeric(1)), c(100,150))
    expect_equal(vapply(values, function(v) v[["SRGF"]], numeric(1)), c(1,1))
  }
})

test_that("socket workers can run formatting helpers without fork", {
  cl <- parallel::makePSOCKcluster(2L); on.exit(parallel::stopCluster(cl))
  parallel::clusterExport(cl, c(".format_wth_value", ".weather_wind_run"), envir=environment())
  actual <- parallel::parLapply(cl, c(13, -99), function(x) .format_wth_value(.weather_wind_run(x)))
  expect_equal(unlist(actual), c("  1123", "   -99"))
})

test_that("native macOS and Windows never select fork workers", {
  expect_false(.soil_use_fork(4L, "unix", "Darwin"))
  expect_false(.soil_use_fork(4L, "windows", "Windows"))
  expect_false(.soil_use_fork(1L, "unix", "Linux"))
  expect_true(.soil_use_fork(4L, "unix", "Linux"))
})

test_that("Alderman conductivity does not shift following fields", {
  layers <- data.frame(SLB=c(20,50), SLMH="-99", SLLL=.1, SDUL=.2, SSAT=.4,
    SRGF=1, SSKS=c("99.996","150"), SBDM=1.4, SLOC=1.2, SLCL=20, SLSI=30,
    SLCF=0, SLNI=-99, SLHW=6, SLHB=-99, SCEC=-99, SADC=-99)
  profile <- list(profile_id="TEST", site="TEST", country="USA", latitude=50, longitude=10,
    scs_family="-99", scom="BN", salb=.13, slu1=6, sldr=.6, slro=73, slnf=1,
    slpf=1, smhb="IB001", smpx="IB001", smke="IB001", layers=layers)
  out <- tempfile(); dir.create(out); on.exit(unlink(out, recursive=TRUE))
  write_dssat_soil_file(profile, out)
  lines <- readLines(file.path(out, "TEST.SOL")); index <- grep("^@  SLB", lines)
  rows <- lines[seq.int(index+1L,length(lines))]; rows <- rows[nzchar(trimws(rows))]
  values <- lapply(rows, function(line) native_values(lines[index], line))
  expect_equal(vapply(values, function(v) v[["SSKS"]], numeric(1)), c(100,150))
  expect_equal(vapply(values, function(v) v[["SBDM"]], numeric(1)), c(1.4,1.4))
})

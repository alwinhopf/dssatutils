.prism_test_records <- function(sdate = "2020-01-01", edate = "2020-12-31") {
  cbind(format(seq(as.Date(sdate), as.Date(edate), by = "day"), "%Y-%m-%d"), "26", "16", "T")
}
.prism_test_response <- function(mat) {
  txt <- jsonlite::toJSON(list(data = mat), auto_unbox = TRUE, digits = NA)
  structure(list(status_code = 200L, content = charToRaw(txt),
                 headers = list(`Content-Type` = "application/json")), class = "response")
}
.prism_test_run <- function(tmp, lat = 30.5, lon = -84.5, pid = "P1", ...) {
  process_weather_prism(data.frame(ID = pid, LAT = lat, LONG = lon), 2020, 2020,
    file.path(tmp, "wth"), "ID", "LAT", "LONG", 1, file.path(tmp, "prism.log"),
    file.path(tmp, "cache"), ...)
}

test_that("ACIS character response coercion uses the shared offline fixture", {
  mat <- jsonlite::fromJSON(test_path("..", "fixtures", "prism_acis_records.json"))$data
  df <- .prism_parse_acis_matrix(mat, 30)
  expect_equal(df$DATE, c("2020001", "2020002", "2020005", "2020006"))
  expect_equal(df$RAIN, c(0, 5.4, -99, -99))
  expect_equal(df$TMAX[2], 27.5)
  expect_true(all(is.finite(df$SRAD) & df$SRAD > 0))
  expect_true(all(.prism_parse_acis_matrix(mat, 30, "none")$SRAD == -99))
})

test_that("ACIS cache reused offline and moved IDs get new weather", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  calls <- 0L
  mat <- .prism_test_records()
  local_mocked_bindings(POST = function(...) {
    calls <<- calls + 1L
    .prism_test_response(mat)
  }, .package = "httr")
  .prism_test_run(tmp)
  expect_equal(calls, 1L)
  expect_true(is_wth_valid(file.path(tmp, "wth/P1.WTH"), start_year = 2020, end_year = 2020,
                          required_columns = c("TMAX", "TMIN", "RAIN", "SRAD")))
  .prism_test_run(tmp, pid = "RENAMED")
  expect_equal(calls, 1L)
  mat[, 2] <- "35"
  .prism_test_run(tmp, lat = 40.5, lon = -100.5)
  expect_equal(calls, 2L)
  expect_true(any(grepl("35.0", readLines(file.path(tmp, "wth/P1.WTH")), fixed = TRUE)))
  expect_length(list.files(file.path(tmp, "cache/acis"), pattern = "\\.json$"), 2)
  expect_length(list.files(file.path(tmp, "cache/acis"), pattern = "^\\.prism-", all.files = TRUE), 0)
})

test_that("ACIS incomplete responses are retried without cache or WTH publication", {
  local_mocked_bindings(Sys.sleep = function(...) NULL, .package = "base")
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  good <- .prism_test_records()
  gap <- good[-21, ]
  duplicate <- good; duplicate[21, ] <- good[20, ]
  reordered <- good; reordered[c(21,22), ] <- good[c(22,21), ]
  missing_temp <- good; missing_temp[21, 2] <- "M"
  missing_rain <- good; missing_rain[21, 4] <- "-999"
  nonfinite <- good; nonfinite[21, 2] <- "Inf"
  malformed <- good[, 1, drop = FALSE]
  cases <- list(short = good[1:2, ], gap = gap, duplicate = duplicate,
                reordered = reordered, missing_temp = missing_temp, missing_rain = missing_rain,
                nonfinite = nonfinite, malformed = malformed, empty = list())
  calls <- 0L
  raw <- NULL
  local_mocked_bindings(POST = function(...) {
    calls <<- calls + 1L
    .prism_test_response(raw)
  }, .package = "httr")
  for (name in names(cases)) {
    raw <- cases[[name]]; calls <- 0L
    dest <- file.path(tmp, name); dir.create(dest)
    .prism_test_run(dest)
    expect_equal(calls, 3L, info = name)
    expect_false(file.exists(file.path(dest, "wth/P1.WTH")), info = name)
    expect_length(list.files(file.path(dest, "cache"), pattern = "\\.json$", recursive = TRUE), 0)
    expect_true(any(grepl("incomplete or invalid ACIS history", readLines(file.path(dest, "prism.log")))))
  }
})

test_that("ACIS legacy and mismatched caches are retried and replaced", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  good <- .prism_test_records("2020-01-01", "2020-01-03")
  request <- .prism_acis_request(30.5, -84.5, "2020-01-01", "2020-01-03")
  wrong_coords <- request; wrong_coords$loc <- "-100.5000,40.5000"
  wrong_units <- request; wrong_units$units <- "degreeF,degreeF,inch"
  cases <- list(legacy = list(data = good), coordinates = list(request = wrong_coords, data = good),
                units = list(request = wrong_units, data = good),
                gap = list(request = request, data = good[-2, ]),
                empty = list(request = request, data = list()), corrupt = "{broken", scalar = '"invalid"')
  calls <- 0L
  local_mocked_bindings(POST = function(...) {
    calls <<- calls + 1L
    .prism_test_response(good)
  }, .package = "httr")
  for (name in names(cases)) {
    f <- file.path(tmp, paste0(name, ".json"))
    writeLines(if (name %in% c("corrupt", "scalar")) cases[[name]] else
                 jsonlite::toJSON(cases[[name]], auto_unbox = TRUE), f)
    calls <- 0L
    result <- .prism_fetch_acis_point(30.5, -84.5, "2020-01-01", "2020-01-03", f)
    expect_equal(result, good, info = name)
    expect_equal(calls, 1L, info = name)
    expect_identical(jsonlite::fromJSON(f)$request, request)
    # Second read must not access HTTP.
    result <- .prism_fetch_acis_point(30.5, -84.5, "2020-01-01", "2020-01-03", f)
    expect_equal(result, good)
    expect_equal(calls, 1L)
  }
})

test_that("ACIS partial response retries before cache publication", {
  local_mocked_bindings(Sys.sleep = function(...) NULL, .package = "base")
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  good <- .prism_test_records("2020-01-01", "2020-01-03")
  calls <- 0L
  local_mocked_bindings(POST = function(...) {
    calls <<- calls + 1L
    .prism_test_response(if (calls == 1L) good[1, , drop = FALSE] else good)
  }, .package = "httr")
  f <- file.path(tmp, "cache.json")
  expect_equal(.prism_fetch_acis_point(30.5, -84.5, "2020-01-01", "2020-01-03", f), good)
  expect_equal(calls, 2L)
  expect_equal(jsonlite::fromJSON(f)$data, good)
})

test_that("ACIS failures preserve existing weather files", {
  local_mocked_bindings(Sys.sleep = function(...) NULL, .package = "base")
  local_mocked_bindings(POST = function(...) .prism_test_response(list()), .package = "httr")
  tmp <- tempfile(); dir.create(tmp); dir.create(file.path(tmp, "wth"))
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  f <- file.path(tmp, "wth/P1.WTH")
  writeLines("previous successful weather", f)
  expect_equal(.prism_test_run(tmp), 0)
  expect_equal(readLines(f), "previous successful weather")
})

test_that("ACIS assembler independently rejects incomplete history", {
  local_mocked_bindings(.prism_fetch_acis_point = function(...) .prism_test_records()[1:2, ],
                       .package = "dssatutils")
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_equal(.prism_test_run(tmp), 0)
  expect_false(file.exists(file.path(tmp, "wth/P1.WTH")))
})

test_that("ACIS none radiation permits only optional missing forcing", {
  local_mocked_bindings(POST = function(...) .prism_test_response(.prism_test_records()), .package = "httr")
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  expect_equal(.prism_test_run(tmp, srad_method = "none"), 1)
  f <- file.path(tmp, "wth/P1.WTH")
  expect_true(is_wth_valid(f, start_year = 2020, end_year = 2020, required_columns = c("TMAX", "TMIN", "RAIN")))
  expect_false(is_wth_valid(f, required_columns = "SRAD"))
})

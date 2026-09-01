library(testthat)

test_that("GRIDMET retains invalid, duplicate and out-of-order points", {
  points <- read.csv(test_path("..", "fixtures", "gridmet_point_alignment.csv"))
  r <- terra::rast(nrows = 1, ncols = 3, nlyrs = 2)
  terra::values(r) <- cbind(c(10, 20, 30), c(11, 21, 31))
  result <- dssatutils:::.gridmet_extract_cells(r, points$cell)
  expect_equal(result, unname(as.matrix(points[c("day1", "day2")])) )
  expect_true(all(is.na(dssatutils:::.gridmet_extract_cells(r, NA_integer_))))
  expect_equal(dssatutils:::.gridmet_extract_cells(r[[1]], 2), matrix(20, 1, 1))
})

test_that("soil rebuild preserves columns and rejects historic shifted caches", {
  fixture <- test_path("..", "fixtures", "soil_mapping_rebuild.csv")
  for (source in c("SSURGO", "GNATSGO")) {
    destination <- tempfile()
    records <- dssatutils::rebuild_soil_files_from_mapping(fixture, destination, source)
    expect_identical(records$ID, "00000001")
    path <- records$path[1]
    expect_null(dssatutils::soil_file_issue(path))
    lines <- readLines(path)
    expect_equal(trimws(substr(lines[4], 2, 11)), "00000001")
    expect_equal(toupper(trimws(substr(lines[4], 14, 24))), source)
    expect_equal(as.numeric(substr(lines[4], 32, 36)), 200)
    expect_equal(as.numeric(substr(lines[6], 26, 33)), 33.682)
    expect_equal(as.numeric(substr(lines[6], 35, 42)), -89.478)
    h <- grep("^@  SLB", lines)
    rows <- lines[h + 1:4]
    expect_equal(as.integer(substr(rows, 1, 6)), c(5L, 20L, 100L, 200L))
    expect_equal(as.numeric(substr(rows, 14, 18)), c(.1, .11, .12, .13))
    expect_equal(as.numeric(substr(rows[1], 38, 42)), 120.5)
    expect_equal(as.numeric(substr(rows[1], 50, 54)), 1.5)
    expect_error(dssatutils::rebuild_soil_files_from_mapping(fixture, destination, source), "new directory")
    expect_identical(readLines(path), lines)
    lines[h + 2:4] <- paste0(" ", lines[h + 2:4])
    writeLines(lines, path)
    expect_match(dssatutils::soil_file_issue(path), "fixed-width")
  }
})

test_that("soil preflight rejects missing files and duplicate depths", {
  expect_identical(dssatutils::soil_file_issue(tempfile()), "SOIL.SOL is missing")
  records <- dssatutils::rebuild_soil_files_from_mapping(
    test_path("..", "fixtures", "soil_mapping_rebuild.csv"), tempfile(), "SSURGO")
  lines <- readLines(records$path[1])
  h <- grep("^@  SLB", lines)
  lines[h + 2] <- lines[h + 1]
  writeLines(lines, records$path[1])
  expect_match(dssatutils::soil_file_issue(records$path[1]), "not strictly increasing")
})

test_that("Alderman retains full-width horizon names, conductivity and sentinels", {
  layers <- read.csv(test_path("..", "fixtures", "soil_mapping_rebuild.csv"))
  names(layers)[match(c("bulk_density", "clay_pct", "silt_pct"), names(layers))] <- c("SBDM", "SLCL", "SLSI")
  layers$SLB <- c(5, 20, 100, 200)
  layers$SLMH <- "Ap/Bt"
  layers$SRGF <- 1
  layers$SSKS <- c(10.08, 120.5, .25, -99)
  layers$SLOC <- layers$om_pct / 1.724
  for (name in c("SLCF", "SLNI", "SLHW", "SLHB", "SCEC", "SADC")) layers[[name]] <- NA_real_
  profile <- list(profile_id="00000001", site="00000001", country="USA", latitude=33.682,
                  longitude=-89.478, scs_family="", scom="SC", salb=.13, slu1=6, sldr=.6,
                  slro=73, slnf=1, slpf=1, smhb="IB001", smpx="IB001", smke="IB001", layers=layers)
  destination <- tempfile()
  dir.create(destination)
  dssatutils:::write_dssat_soil_file(profile, destination)
  path <- file.path(destination, "00000001.SOL")
  expect_null(dssatutils::soil_file_issue(path))
  lines <- readLines(path)
  h <- grep("^@  SLB", lines)
  rows <- lines[h + 1:4]
  expect_equal(substr(rows, 8, 12), rep("Ap/Bt", 4))
  expect_equal(as.numeric(substr(rows, 38, 42)), c(10.08, 120.5, .25, -99))
  lines[h + 1:4] <- substring(lines[h + 1:4], 2)
  writeLines(lines, path)
  expect_match(dssatutils::soil_file_issue(path), "fixed-width")
})

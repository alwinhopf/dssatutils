library(testthat)
library(dssatutils)

# Helper: assert a .WTH file is valid
assert_wth_valid <- function(path) {
  expect_true(file.exists(path))
  content <- readLines(path)
  expect_true(any(grepl("\\$WEATHER DATA", content)))
  expect_true(any(grepl("@  DATE", content)))
}

# Helper: assert a .SOL file is valid
assert_sol_valid <- function(path) {
  expect_true(file.exists(path))
  content <- readLines(path)
  has_header <- any(grepl("\\*SOILS:", content)) || any(grepl("^\\*", content))
  expect_true(has_header)
  expect_true(any(grepl("@  SLB", content)))
}

# ===================================================================
# WEATHER SOURCES
# ===================================================================

test_that("process_weather_openmeteo runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")
  
  local_mocked_bindings(
    GET = function(...) {
      mock_response <- list(
        status_code = 200,
        content = charToRaw('{"daily": {
          "time": ["2010-01-01", "2010-01-02"],
          "temperature_2m_max": [25.0, 26.0],
          "temperature_2m_min": [15.0, 16.0],
          "precipitation_sum": [0.0, 1.2],
          "shortwave_radiation_sum": [18.0, 19.0],
          "wind_speed_10m_max": [3.0, 4.0]
        }}')
      )
      class(mock_response) <- "response"
      mock_response
    },
    status_code = function(r) r$status_code,
    content = function(r, as = "text", encoding = "UTF-8") rawToChar(r$content),
    .package = "httr"
  )
  
  process_weather_openmeteo(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

test_that("process_weather_nasapower runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")
  
  local_mocked_bindings(
    get_power = function(community, lonlat, pars, dates, temporal_api) {
      data.frame(
        YEAR = c(2010, 2010),
        MM = c(1, 1),
        DOY = c(1, 2),
        T2M_MAX = c(25.0, 26.0),
        T2M_MIN = c(15.0, 16.0),
        ALLSKY_SFC_SW_DWN = c(18.0, 19.0),
        PRECTOTCORR = c(0.0, 1.2),
        T2MDEW = c(12.0, 13.0),
        RH2M = c(80.0, 82.0),
        WS2M = c(3.0, 4.0)
      )
    },
    .package = "nasapower"
  )
  
  process_weather_nasapower(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

test_that("process_weather_daymet runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")
  
  local_mocked_bindings(
    download_daymet = function(lat, lon, start, end, internal, silent) {
      list(
        data = data.frame(
          year = c(2010, 2010),
          yday = c(1, 2),
          tmax..deg.c. = c(25.0, 26.0),
          tmin..deg.c. = c(15.0, 16.0),
          prcp..mm.day. = c(0.0, 1.2),
          srad..W.m.2. = c(200.0, 210.0),
          dayl..s. = c(43200, 43200),
          vp..Pa. = c(1200, 1250)
        )
      )
    },
    .package = "daymetr"
  )
  
  process_weather_daymet(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

test_that("process_weather_nasapower_chirps runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")
  
  local_mocked_bindings(
    vect = function(...) "mock_vect",
    project = function(...) "mock_project",
    crs = function(...) "mock_crs",
    extract = function(...) data.frame(lyr1 = 5.0, lyr2 = 6.0),
    rast = function(...) "mock_rast",
    time = function(...) as.Date(c("2010-01-01", "2010-01-02")),
    .package = "terra"
  )
  
  local_mocked_bindings(
    download.file = function(...) 0,
    .package = "utils"
  )
  
  local_mocked_bindings(
    get_power = function(community, lonlat, pars, dates, temporal_api) {
      data.frame(
        YEAR = c(2010, 2010),
        MM = c(1, 1),
        DOY = c(1, 2),
        T2M_MAX = c(25.0, 26.0),
        T2M_MIN = c(15.0, 16.0),
        ALLSKY_SFC_SW_DWN = c(18.0, 19.0),
        PRECTOTCORR = c(0.0, 1.2),
        T2MDEW = c(12.0, 13.0),
        RH2M = c(80.0, 82.0),
        WS2M = c(3.0, 4.0)
      )
    },
    .package = "nasapower"
  )
  
  process_weather_nasapower_chirps(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file,
    chirps_cache_dir = file.path(work_dir, "chirps_cache")
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

test_that("process_weather_agera5 runs successfully with mocks", {
  # AgERA5 is an optional source: ecmwfr is in Suggests, so it may be absent on
  # CI. Mocking requires the namespace to load, hence skip when not installed.
  skip_if_not_installed("ecmwfr")

  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))

  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")

  local_mocked_bindings(
    wf_request = function(...) 0,
    .package = "ecmwfr"
  )
  
  local_mocked_bindings(
    vect = function(...) "mock_vect",
    extract = function(...) matrix(c(298.15, 299.15), nrow = 1),
    rast = function(...) "mock_rast",
    time = function(...) as.Date(c("2010-01-01", "2010-01-02")),
    .package = "terra"
  )
  
  process_weather_agera5(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file,
    agera5_cache_dir = file.path(work_dir, "agera5_cache")
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

test_that("process_weather_gridmet runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")
  
  local_mocked_bindings(
    GET = function(...) {
      mock_response <- list(status_code = 200)
      class(mock_response) <- "response"
      mock_response
    },
    write_disk = function(...) "mock_disk",
    .package = "httr"
  )
  
  local_mocked_bindings(
    rast = function(...) {
      structure(list(), class = "MockSpatRaster")
    },
    cellFromXY = function(...) 1,
    nlyr = function(...) 365,
    ext = function(...) c(-125, -66.5, 25, 49.5),
    extract = function(...) matrix(c(280.0, 281.0), nrow = 1),
    .package = "terra"
  )
  
  `[.MockSpatRaster` <- function(x, i, j, ...) {
    matrix(c(280.0, 281.0), nrow = 1)
  }
  registerS3method("[", "MockSpatRaster", `[.MockSpatRaster`, envir = parent.frame())
  
  local_mocked_bindings(
    st_transform = function(x, ...) x,
    st_as_sf = function(x, ...) x,
    st_coordinates = function(x, ...) matrix(c(-90.0, 40.0), nrow = 1),
    .package = "sf"
  )
  
  # GridMET downloads files into gridmet_cache_dir and checks if they exist.
  # Let's create dummy files for all 4 variables so that the existence check passes and downloading is bypassed/mocked.
  cache_dir <- file.path(work_dir, "gridmet_cache")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  for (var in c("tmmn", "tmmx", "pr", "srad")) {
    writeLines("", file.path(cache_dir, paste0(var, "_2010.nc")))
  }
  
  process_weather_gridmet(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file,
    gridmet_cache_dir = cache_dir
  )
  
  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

# ===================================================================
# SOIL SOURCES
# ===================================================================

test_that("process_soils_soilgrids runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  # Create a dummy master SOL file
  master_sol <- file.path(work_dir, "master.SOL")
  writeLines(c(
    "*SOILS: Dummy Master",
    "*TEST0001     TEST_SOIL       40.000   -90.000",
    "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
    " TEST_SOIL   World         40.000   -90.000",
    "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
    "    BN   .13     6    .6    73     1     1 IB001 IB001 IB001",
    "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC",
    "     5   -99 0.100 0.200 0.300  1.00  10.0  1.40  1.00  20.0  40.0   0.0   -99   -99   -99   -99   -99"
  ), master_sol)
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  output_csv <- file.path(work_dir, "soil_map.csv")
  output_sol_dir <- file.path(work_dir, "individual_sol")
  
  local_mocked_bindings(
    st_as_sf = function(x, coords, crs, ...) {
      class(x) <- c("sf", "data.frame")
      x
    },
    st_transform = function(x, ...) x,
    st_coordinates = function(x, ...) {
      matrix(c(-90.0, 40.0), nrow = nrow(x))
    },
    st_nearest_feature = function(x, y) {
      rep(1, nrow(x))
    },
    .package = "sf"
  )
  
  process_soils_soilgrids(
    grid_points = shapefile,
    source_sol_file = master_sol,
    output_csv_path = output_csv,
    output_sol_dir = output_sol_dir,
    id_col = "ID",
    numeric_only_ids = FALSE
  )
  
  expect_true(file.exists(output_csv))
  assert_sol_valid(file.path(output_sol_dir, "TEST1.SOL"))
})

test_that("process_soils_soilgrids_online runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  # Inject USE_REST_API = TRUE in the global scope for the test
  assign("USE_REST_API", TRUE, envir = globalenv())
  on.exit(rm("USE_REST_API", envir = globalenv()), add = TRUE)
  
  # Mock JSON REST response
  local_mocked_bindings(
    GET = function(...) {
      mock_response <- list(
        status_code = 200,
        content = charToRaw('{"properties": {"layers": [
          {"name": "clay", "depths": [{"label": "0-5cm", "values": {"mean": 200}}]},
          {"name": "sand", "depths": [{"label": "0-5cm", "values": {"mean": 400}}]},
          {"name": "silt", "depths": [{"label": "0-5cm", "values": {"mean": 400}}]},
          {"name": "soc", "depths": [{"label": "0-5cm", "values": {"mean": 150}}]},
          {"name": "bdod", "depths": [{"label": "0-5cm", "values": {"mean": 130}}]},
          {"name": "cfvo", "depths": [{"label": "0-5cm", "values": {"mean": 10}}]}
        ]}}')
      )
      class(mock_response) <- "response"
      mock_response
    },
    status_code = function(r) r$status_code,
    content = function(r, as = "text", encoding = "UTF-8") rawToChar(r$content),
    .package = "httr"
  )
  
  # sf methods mock
  local_mocked_bindings(
    st_transform = function(x, ...) x,
    st_coordinates = function(x, ...) matrix(c(-90.0, 40.0), nrow = 1),
    .package = "sf"
  )
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  shapefile[["geometry"]] <- sf::st_sfc(sf::st_point(c(-90.0, 40.0)))
  class(shapefile) <- c("sf", "data.frame")
  attr(shapefile, "sf_column") <- "geometry"
  
  output_csv <- file.path(work_dir, "soil_map.csv")
  output_sol_dir <- file.path(work_dir, "individual_sol")
  dir.create(output_sol_dir, recursive = TRUE, showWarnings = FALSE)
  
  process_soils_soilgrids_online(
    gridfile = shapefile,
    soilfile_csv_path = output_csv,
    output_sol_dir = output_sol_dir,
    id_col = "ID"
  )
  
  expect_true(file.exists(output_csv))
  assert_sol_valid(file.path(output_sol_dir, "TEST1.SOL"))
})

test_that("process_soils_ssurgo runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  shapefile[["geometry"]] <- sf::st_sfc(sf::st_point(c(-90.0, 40.0)))
  class(shapefile) <- c("sf", "data.frame")
  attr(shapefile, "sf_column") <- "geometry"
  
  local_mocked_bindings(
    robust_SDA_spatialQuery = function(point_sf, what, ...) {
      list(ok = TRUE, data = data.frame(mukey = "12345"), error = NA_character_)
    },
    robust_SDA_query = function(query, ...) {
      if (grepl("brockdepmin", query)) {
        list(ok = TRUE, data = data.frame(mukey = "12345", brockdepmin = 200.0), error = NA_character_)
      } else {
        list(ok = TRUE, data = data.frame(
          mukey = "12345",
          cokey = "cokey1",
          comppct_r = 100,
          hzdept_r = 0,
          hzdepb_r = 15,
          claytotal_r = 20.0,
          sandtotal_r = 40.0,
          om_r = 1.5,
          dbthirdbar_r = 1.4
        ), error = NA_character_)
      }
    },
    .package = "dssatutils"
  )
  
  local_mocked_bindings(
    st_as_sf = function(...) {
      sf_obj <- data.frame(ID = "TEST1")
      sf_obj[["geometry"]] <- sf::st_sfc(sf::st_point(c(-90.0, 40.0)))
      class(sf_obj) <- c("sf", "data.frame")
      attr(sf_obj, "sf_column") <- "geometry"
      sf_obj
    },
    st_drop_geometry = function(x) {
      x[["geometry"]] <- NULL
      class(x) <- "data.frame"
      x
    },
    st_coordinates = function(...) matrix(c(-90.0, 40.0), nrow = 1),
    .package = "sf"
  )
  
  output_csv <- file.path(work_dir, "soil_map.csv")
  output_sol_dir <- file.path(work_dir, "individual_sol")
  dir.create(output_sol_dir, recursive = TRUE, showWarnings = FALSE)
  
  process_soils_ssurgo(
    grid_points = shapefile,
    output_dir_csv = output_csv,
    output_dir_individual = output_sol_dir,
    n_cores = 1,
    id_col = "ID",
    lat_col = "LAT",
    long_col = "LONG",
    format_sql_func = function(x) paste0("('", x, "')")
  )
  
  expect_true(file.exists(output_csv))
  assert_sol_valid(file.path(output_sol_dir, "TEST1.SOL"))
})

test_that("process_soils_ssurgo_alderman runs successfully with mocks", {
  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  shapefile[["geometry"]] <- sf::st_sfc(sf::st_point(c(-90.0, 40.0)))
  class(shapefile) <- c("sf", "data.frame")
  attr(shapefile, "sf_column") <- "geometry"
  
  local_mocked_bindings(
    robust_SDA_spatialQuery = function(point_sf, what, ...) {
      data.frame(mukey = "12345")
    },
    robust_SDA_query = function(query, ...) {
      if (grepl("muaggatt", query)) {
        data.frame(mukey = "12345", brockdepmin = 200.0)
      } else if (grepl("component", query)) {
        data.frame(
          compname = "Miami",
          cokey = "cokey1",
          mukey = "12345",
          comppct_r = 100,
          hydgrp = "B",
          slope_r = 2.0,
          drainage = "Well drained",
          albedodry_r = 0.13
        )
      } else {
        data.frame(
          hzdept_r = 0,
          hzdepb_r = 15,
          dbovendry_r = 1.45,
          dbtenthbar_r = 1.4,
          dbthirdbar_r = 1.4,
          dbfifteenbar_r = 1.5,
          wsatiated_r = 45.0,
          wtenthbar_r = 25.0,
          wthirdbar_r = 20.0,
          partdensity = 2.65,
          ksat_r = 15.0,
          wfifteenbar_r = 10.0,
          sandtotal_r = 40.0,
          claytotal_r = 20.0,
          silttotal_r = 40.0,
          om_r = 1.5,
          hzname = "Ap",
          fragvol_r = NA_real_,
          cokey = "cokey1"
        )
      }
    },
    .package = "dssatutils"
  )
  
  local_mocked_bindings(
    st_as_sf = function(...) {
      sf_obj <- data.frame(ID = "TEST1")
      sf_obj[["geometry"]] <- sf::st_sfc(sf::st_point(c(-90.0, 40.0)))
      class(sf_obj) <- c("sf", "data.frame")
      attr(sf_obj, "sf_column") <- "geometry"
      sf_obj
    },
    st_drop_geometry = function(x) {
      x[["geometry"]] <- NULL
      class(x) <- "data.frame"
      x
    },
    st_coordinates = function(...) matrix(c(-90.0, 40.0), nrow = 1),
    st_as_text = function(...) "POINT(-90 40)",
    .package = "sf"
  )
  
  output_csv <- file.path(work_dir, "soil_map.csv")
  output_sol_dir <- file.path(work_dir, "individual_sol")
  dir.create(output_sol_dir, recursive = TRUE, showWarnings = FALSE)
  
  process_soils_ssurgo_alderman(
    grid_points = shapefile,
    output_dir_csv = output_csv,
    output_dir_individual = output_sol_dir,
    n_cores = 1,
    id_col = "ID",
    lat_col = "LAT",
    long_col = "LONG",
    format_sql_func = function(x) paste0("('", x, "')")
  )
  
  expect_true(file.exists(output_csv))
  assert_sol_valid(file.path(output_sol_dir, "TEST1.SOL"))
})

test_that("process_soils_hwsd runs successfully with mocks", {
  # HWSD uses DBI + RSQLite, both in Suggests, so they may be absent on CI.
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")

  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))

  # Create a dummy raster TIFF file
  dummy_tif <- file.path(work_dir, "dummy_hwsd.tif")
  writeLines("", dummy_tif)
  
  # Create a dummy SQLite DB
  dummy_db <- file.path(work_dir, "dummy_hwsd.sqlite")
  con <- DBI::dbConnect(RSQLite::SQLite(), dummy_db)
  DBI::dbExecute(con, "CREATE TABLE HWSD2_LAYERS (
    HWSD2_SMU_ID INTEGER, SEQUENCE INTEGER, SHARE REAL,
    TOPDEP REAL, BOTDEP REAL, SAND REAL, CLAY REAL, SILT REAL,
    BULK REAL, ORG_CARBON REAL, COARSE REAL
  )")
  DBI::dbExecute(con, "INSERT INTO HWSD2_LAYERS VALUES (1, 1, 100, 0, 15, 40.0, 20.0, 40.0, 1.40, 1.0, 0.0)")
  DBI::dbDisconnect(con)
  
  local_mocked_bindings(
    rast = function(...) "mock_rast",
    vect = function(...) "mock_vect",
    same.crs = function(...) TRUE,
    extract = function(...) data.frame(lyr1 = 1),
    .package = "terra"
  )
  
  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  output_csv <- file.path(work_dir, "soil_map.csv")
  output_sol_dir <- file.path(work_dir, "individual_sol")
  
  process_soils_hwsd(
    grid_points = shapefile,
    hwsd_raster_file = dummy_tif,
    hwsd_db_file = dummy_db,
    output_csv_path = output_csv,
    output_sol_dir = output_sol_dir,
    id_col = "ID",
    lat_col = "LAT",
    long_col = "LONG"
  )
  
  expect_true(file.exists(output_csv))
  assert_sol_valid(file.path(output_sol_dir, "TEST1.SOL"))
})

test_that("process_weather_era5_land runs successfully with mocks", {
  skip_if_not_installed("ecmwfr")

  work_dir <- tempfile()
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE))

  shapefile <- data.frame(ID = "TEST1", LAT = 40.0, LONG = -90.0)
  log_file <- file.path(work_dir, "error.log")

  # Mock downloading by writing a synthetic CSV to the raw_csv destination
  local_mocked_bindings(
    .download_era5_land_point_csv = function(latitude, longitude, start_date, end_date, target_file, ...) {
      df <- data.frame(
        time = seq(as.POSIXct("2010-01-01 00:00:00", tz = "UTC"), as.POSIXct("2010-01-02 23:00:00", tz = "UTC"), by = "hour"),
        `2m_temperature` = runif(48, 280, 295),
        `2m_dewpoint_temperature` = runif(48, 275, 285),
        `total_precipitation` = runif(48, 0, 0.005),
        `surface_solar_radiation_downwards` = runif(48, 1e5, 1e6),
        `10m_u_component_of_wind` = runif(48, -5, 5),
        `10m_v_component_of_wind` = runif(48, -5, 5),
        check.names = FALSE
      )
      readr::write_csv(df, target_file)
    },
    .package = "dssatutils"
  )

  process_weather_era5_land(
    shapefile = shapefile,
    start_year = 2010,
    end_year = 2010,
    output_dir = work_dir,
    id_col = "ID",
    lat_col = "LAT",
    lon_col = "LONG",
    n_cores = 1,
    log_file = log_file
  )

  assert_wth_valid(file.path(work_dir, "TEST1.WTH"))
})

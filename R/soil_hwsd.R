# File: soil_hwsd.R
# ---------------------------------------------------------------------------
# Soil source: FAO/IIASA Harmonized World Soil Database v2.0 (HWSD2).
#
# HWSD2 is the FAO "official" harmonized global soil database (~1 km), the
# long-standing reference for GLOBAL gridded crop-model studies. It is NOT a
# streaming API: download it once from FAO and point the pipeline at two files
# (mirrors the SOILGRIDS_10K external-file model):
#   * hwsd_raster_file : HWSD2 raster of mapping-unit (SMU) IDs (GeoTIFF/BIL)
#   * hwsd_db_file     : HWSD2 attribute database (SQLite)
#   FAO HWSD v2.0: https://www.fao.org/soils-portal/data-hub/soil-maps-and-databases/harmonized-world-soil-database-v2-0/
#
# Samples the raster at each point to get the SMU ID, looks up the dominant
# soil component's layers in the SQLite DB, computes DSSAT physics (Saxton &
# Rawls 2006, shared with the SoilGrids module), and writes per-point .SOL
# files + a mapping CSV. Points over no-data cells are skipped with a warning.
#
# NOTE: validated structurally; run once against your FAO HWSD2 download to
# confirm the SQLite column names match (handled tolerantly below).
# Requires: terra, RSQLite, DBI, dplyr. Reuses calculate_soil_physics() and
# format_dssat_sol_file() from soil_soilgrids_online.R (source that first).
# ---------------------------------------------------------------------------

library(terra)
library(dplyr)
# DBI + RSQLite are loaded lazily inside process_soils_hwsd() so that simply
# sourcing this file (which the main pipeline does for ALL soil modules at
# startup) does not require the SQLite packages unless you actually run HWSD.

# Tolerant column matching for the HWSD2 layer table.
.hwsd_pick <- function(cols, candidates) {
  lc <- tolower(cols)
  for (cand in candidates) {
    hit <- which(lc == tolower(cand))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  NA_character_
}

process_soils_hwsd <- function(grid_points, hwsd_raster_file, hwsd_db_file,
                               output_csv_path, output_sol_dir,
                               id_col = "ID", lat_col = "LAT", long_col = "LONG") {
  message("--- Starting HWSD2 Extraction ---")
  # Load the SQLite packages only now (HWSD-only deps; see note at top of file).
  for (pkg in c("DBI", "RSQLite")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("HWSD needs the '%s' package. install.packages('%s')", pkg, pkg))
  }
  for (f in c(hwsd_raster_file, hwsd_db_file)) {
    if (!file.exists(f)) stop(sprintf("HWSD2 file not found: %s", f))
  }
  if (!dir.exists(output_sol_dir)) dir.create(output_sol_dir, recursive = TRUE)

  ids  <- as.character(grid_points[[id_col]])
  # Coordinates: prefer sf/terra geometry, else explicit lat/long columns.
  if (inherits(grid_points, "sf")) {
    g <- sf::st_transform(grid_points, 4326)
    co <- sf::st_coordinates(g); lons <- co[, 1]; lats <- co[, 2]
  } else {
    lons <- as.numeric(grid_points[[long_col]]); lats <- as.numeric(grid_points[[lat_col]])
  }

  # 1. Sample HWSD2 SMU IDs at each point.
  r <- terra::rast(hwsd_raster_file)
  pts <- terra::vect(data.frame(lon = lons, lat = lats),
                     geom = c("lon", "lat"), crs = "EPSG:4326")
  if (!terra::same.crs(pts, r)) pts <- terra::project(pts, terra::crs(r))
  smu <- terra::extract(r, pts, ID = FALSE)[, 1]
  smu[is.nan(smu)] <- NA

  # 2. Load the HWSD2 layer table and resolve columns.
  con <- DBI::dbConnect(RSQLite::SQLite(), hwsd_db_file)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tbls <- DBI::dbListTables(con)
  layer_tbl <- .hwsd_pick(tbls, c("HWSD2_LAYERS", "HWSD2_LAYER", "LAYERS", "D_LAYERS"))
  if (is.na(layer_tbl)) stop(sprintf("No HWSD2 layer table found. Tables: %s",
                                     paste(tbls, collapse = ", ")))
  layers <- DBI::dbReadTable(con, layer_tbl)
  cn <- names(layers)
  col <- list(
    smu = .hwsd_pick(cn, c("HWSD2_SMU_ID","SMU_ID","HWSD2_SMU","SMU")),
    seq = .hwsd_pick(cn, c("SEQUENCE","SEQ")),
    share = .hwsd_pick(cn, c("SHARE","PERCENT","PCT")),
    top = .hwsd_pick(cn, c("TOPDEP","TOP_DEPTH","TOP")),
    bot = .hwsd_pick(cn, c("BOTDEP","BOT_DEPTH","BOTTOM","BOT")),
    sand = .hwsd_pick(cn, c("SAND","SAND_PCT")),
    silt = .hwsd_pick(cn, c("SILT","SILT_PCT")),
    clay = .hwsd_pick(cn, c("CLAY","CLAY_PCT")),
    bulk = .hwsd_pick(cn, c("BULK","BULK_DENSITY","BD","REF_BULK_DENSITY")),
    oc = .hwsd_pick(cn, c("ORG_CARBON","OC","ORGANIC_CARBON","SOC")),
    coarse = .hwsd_pick(cn, c("COARSE","GRAVEL","CFVO")))
  if (any(is.na(c(col$smu, col$sand, col$clay, col$bot))))
    stop(sprintf("HWSD2 table '%s' missing required columns. Found: %s",
                 layer_tbl, paste(cn, collapse = ", ")))

  # 3. Build per-point layers (dominant component) + physics + write .SOL.
  success <- 0; errors <- 0; skipped <- character(0)
  log_path <- file.path(output_sol_dir, "soil_processing_errors.log")
  cat(sprintf("Log started: %s\n", Sys.time()), file = log_path)

  for (i in seq_along(ids)) {
    pid <- ids[i]; s <- smu[i]
    if (is.na(s)) { skipped <- c(skipped, pid); next }
    sub <- layers[layers[[col$smu]] == s, , drop = FALSE]
    if (nrow(sub) == 0) { skipped <- c(skipped, pid); next }
    # Dominant component: highest SHARE, else lowest SEQUENCE.
    if (!is.na(col$share) && !is.na(col$seq)) {
      keep_seq <- sub[[col$seq]][which.max(sub[[col$share]])]
      sub <- sub[sub[[col$seq]] == keep_seq, , drop = FALSE]
    } else if (!is.na(col$seq)) {
      sub <- sub[sub[[col$seq]] == min(sub[[col$seq]]), , drop = FALSE]
    }

    site <- lapply(seq_len(nrow(sub)), function(k) {
      lyr <- sub[k, ]
      bot <- as.numeric(lyr[[col$bot]])
      top <- if (!is.na(col$top)) as.numeric(lyr[[col$top]]) else max(0, bot - 20)
      sand <- as.numeric(lyr[[col$sand]]); clay <- as.numeric(lyr[[col$clay]])
      silt <- if (!is.na(col$silt)) as.numeric(lyr[[col$silt]]) else max(0, 100 - sand - clay)
      oc <- if (!is.na(col$oc)) as.numeric(lyr[[col$oc]]) else NA
      phys <- calculate_soil_physics(sand, clay, ifelse(is.na(oc), 1.0, oc * 1.724))
      data.frame(ID = pid, latitude = lats[i], longitude = lons[i],
                 depth_bottom = bot, depth_center = (top + bot) / 2,
                 sand = sand, clay = clay, silt = silt,
                 bdod = if (!is.na(col$bulk)) as.numeric(lyr[[col$bulk]]) else NA,
                 soc_pct = oc,
                 cfvo = if (!is.na(col$coarse)) as.numeric(lyr[[col$coarse]]) else 0,
                 SLLL = phys$SLLL, SDUL = phys$SDUL, SSAT = phys$SSAT)
    })
    site_df <- do.call(rbind, site)

    tryCatch({
      format_dssat_sol_file(site_df, output_sol_dir,
                            source_name = "FAO HWSD v2.0", source_tag = "HWSD2")
      success <- success + 1
    }, error = function(e) {
      errors <<- errors + 1
      cat(sprintf("ID: %s | Error: %s\n", pid, conditionMessage(e)),
          file = log_path, append = TRUE)
    })
  }

  if (length(skipped) > 0)
    message(sprintf("  [warn] %d point(s) had no HWSD2 mapping unit and were skipped: %s%s",
                    length(skipped), paste(head(skipped, 10), collapse = ", "),
                    if (length(skipped) > 10) " ..." else ""))

  # 4. Mapping CSV (ID -> SOIL_ID == ID).
  utils::write.csv(data.frame(ID = ids, SOIL_ID = ids),
                   output_csv_path, row.names = FALSE)
  message(sprintf("HWSD2 processing complete. Success: %d, Errors: %d, Skipped: %d",
                  success, errors, length(skipped)))
}

# ==============================================================================
#  SOIL HELPER: SOILGRIDS / EXTERNAL MASTER .SOL -> per-point .SOL + mapping CSV
#  Filename: soil_soilgrids.R  (FIXED: no top-level execution; numeric-only IDs)
# ==============================================================================
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(stringr)
  library(readr)
})

#' Process SoilGrids / external DSSAT soil master file into per-point SOLs + mapping CSV
#'
#' @param grid_points sf POINT layer (recommended) or data.frame with lon/lat columns.
#' @param source_sol_file Path to master DSSAT .SOL file containing many profiles.
#' @param output_csv_path Path to write mapping CSV (must include SOIL_ID column for pipeline).
#' @param output_sol_dir Directory to write individual <SOIL_ID>.SOL files.
#' @param id_col Column name in grid_points that contains the grid-point ID.
#' @param numeric_only_ids If TRUE, SOIL_ID will be digits only (recommended for DSSAT spatial runs).
#' @param numeric_width Width for zero-padded numeric IDs (default 8 -> 00000028).
#' @return data.frame mapping with columns: ID (point), SOIL_ID, SOURCE_SOIL_ID, lat, lon
process_soils_soilgrids <- function(
  grid_points,
  source_sol_file,
  output_csv_path,
  output_sol_dir,
  id_col = "ID",
  numeric_only_ids = TRUE,
  numeric_width = 8
) {

  message(paste("Parsing external soil file:", source_sol_file))

  if (!file.exists(source_sol_file)) {
    stop(paste("CRITICAL ERROR: The external soil file was not found at:", source_sol_file))
  }

  if (!dir.exists(output_sol_dir)) dir.create(output_sol_dir, recursive = TRUE, showWarnings = FALSE)

  # --------------------------------------------------------------------
  # Helpers
  # --------------------------------------------------------------------
  make_numeric_id <- function(x, width = 8) {
    x <- as.character(x)
    d <- gsub("\\D+", "", x)       # keep digits only
    ifelse(
      nchar(d) == 0, NA_character_,
      ifelse(
        nchar(d) >= width,
        substr(d, nchar(d) - width + 1, nchar(d)),
        paste0(strrep("0", width - nchar(d)), d)
      )
    )
  }

  rewrite_header_id <- function(header_line, new_id) {
    # Replace the first token after '*' with new_id (preserves spacing/format)
    new_id <- as.character(new_id)
    if (nchar(new_id) > 10) new_id <- substr(new_id, nchar(new_id) - 10 + 1, nchar(new_id))
    sub("^\\*\\S+", paste0("*", new_id), header_line)
  }

  # --------------------------------------------------------------------
  # 1. Read the Master File (latin1 to avoid encoding crashes)
  # --------------------------------------------------------------------
  lines <- readLines(source_sol_file, warn = FALSE, encoding = "latin1")

  starts <- grep("^\\*.*", lines)  # profile header lines start with *
  if (length(starts) == 0) stop("No profiles found (no lines starting with '*').")

  ends <- c(tail(starts, -1) - 1, length(lines))

  # Extract source IDs: first token after '*'
  raw_ids <- lines[starts]
  source_ids <- str_extract(raw_ids, "(?<=\\*)[^\\s]+")

  # --------------------------------------------------------------------
  # 2. Extract Coordinates (scan a small chunk for @SITE and lat/lon)
  # --------------------------------------------------------------------
  message("Extracting coordinates (scanning for valid Lat/Lon patterns)...")

  parse_lat_lon_from_chunk <- function(txt_chunk) {
    site_tag_idx <- grep("^@SITE", txt_chunk)
    if (length(site_tag_idx) == 0) return(c(NA_real_, NA_real_))

    candidate_lines <- txt_chunk[(site_tag_idx + 1) : min((site_tag_idx + 3), length(txt_chunk))]

    for (ln in candidate_lines) {
      if (!str_detect(ln, "[0-9]")) next

      raw_matches <- str_extract_all(
        ln,
        "-?\\d+\\.\\d+|-?\\d+\\.|-?\\.\\d+|-?\\d+"
      )[[1]]

      nums <- suppressWarnings(as.numeric(raw_matches))
      candidates <- nums[nums != -99]

      if (length(candidates) >= 2) {
        for (j in 1:(length(candidates) - 1)) {
          lat_c <- candidates[j]
          lon_c <- candidates[j + 1]
          if (!is.na(lat_c) && !is.na(lon_c) && abs(lat_c) <= 90 && abs(lon_c) <= 180) {
            return(c(lat_c, lon_c))
          }
        }
      }
    }
    return(c(NA_real_, NA_real_))
  }

  soil_df_list <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    s <- starts[i]
    e <- ends[i]
    chunk_end <- min(e, s + 25)
    chunk <- lines[s:chunk_end]
    coords <- parse_lat_lon_from_chunk(chunk)

    soil_df_list[[i]] <- data.frame(
      source_soil_id = source_ids[i],
      lat = coords[1],
      lon = coords[2],
      start_line = s,
      end_line = e,
      stringsAsFactors = FALSE
    )
  }

  soil_lib_df <- bind_rows(soil_df_list) %>%
    filter(!is.na(lat) & !is.na(lon))

  if (nrow(soil_lib_df) == 0) {
    stop("No soil profiles with valid lat/lon were found in the master .SOL (check @SITE blocks).")
  }

  soil_sf <- st_as_sf(soil_lib_df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)

  # --------------------------------------------------------------------
  # 3. Prepare grid points as sf (WGS84)
  # --------------------------------------------------------------------
  if (inherits(grid_points, "sf")) {
    gp_sf <- st_transform(grid_points, 4326)
    gp_coords <- st_coordinates(gp_sf)
    gp_df <- gp_sf %>%
      mutate(.gp_lon = gp_coords[,1], .gp_lat = gp_coords[,2])
  } else {
    gp_df <- as.data.frame(grid_points)
    # try common lon/lat names
    lon_col <- intersect(names(gp_df), c("lon","Lon","LON","x","X","LONG","longitude","Longitude"))
    lat_col <- intersect(names(gp_df), c("lat","Lat","LAT","y","Y","LATI","latitude","Latitude"))
    if (length(lon_col) == 0 || length(lat_col) == 0) {
      stop("grid_points must be sf, or a data.frame with lon/lat columns.")
    }
    gp_df$.gp_lon <- gp_df[[lon_col[1]]]
    gp_df$.gp_lat <- gp_df[[lat_col[1]]]
    gp_sf <- st_as_sf(gp_df, coords = c(".gp_lon",".gp_lat"), crs = 4326, remove = FALSE)
  }

  if (!(id_col %in% names(gp_df))) {
    stop(paste("id_col", shQuote(id_col), "not found in grid_points."))
  }

  # --------------------------------------------------------------------
  # 4. Nearest soil profile to each grid point
  # --------------------------------------------------------------------
  nearest_idx <- st_nearest_feature(gp_sf, soil_sf)
  nearest <- soil_lib_df[nearest_idx, , drop = FALSE]

  point_ids_raw <- as.character(gp_df[[id_col]])

  soil_ids <- if (numeric_only_ids) make_numeric_id(point_ids_raw, width = numeric_width) else point_ids_raw
  # fallback if any NA
  na_idx <- which(is.na(soil_ids) | soil_ids == "")
  if (length(na_idx) > 0) {
    soil_ids[na_idx] <- sprintf(paste0("%0", numeric_width, "d"), na_idx)
  }

  mapping <- data.frame(
    ID = point_ids_raw,
    SOIL_ID = soil_ids,
    SOURCE_SOIL_ID = nearest$source_soil_id,
    SOIL_LAT = nearest$lat,
    SOIL_LON = nearest$lon,
    START_LINE = nearest$start_line,
    END_LINE = nearest$end_line,
    stringsAsFactors = FALSE
  )

  # --------------------------------------------------------------------
  # 5. Write mapping CSV
  # --------------------------------------------------------------------
  write_csv(mapping %>% select(ID, SOIL_ID, SOURCE_SOIL_ID, SOIL_LAT, SOIL_LON), output_csv_path)

  # --------------------------------------------------------------------
  # 6. Write per-point .SOL files with rewritten profile IDs
  # --------------------------------------------------------------------
  message("Writing individual .SOL files...")
  for (i in seq_len(nrow(mapping))) {
    s <- mapping$START_LINE[i]
    e <- mapping$END_LINE[i]
    prof_lines <- lines[s:e]

    # rewrite header line
    prof_lines[1] <- rewrite_header_id(prof_lines[1], mapping$SOIL_ID[i])

    out_file <- file.path(output_sol_dir, paste0(mapping$SOIL_ID[i], ".SOL"))
    writeLines(prof_lines, out_file, useBytes = TRUE)
  }

  invisible(mapping %>% select(ID, SOIL_ID, SOURCE_SOIL_ID, SOIL_LAT, SOIL_LON))
}

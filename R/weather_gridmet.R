# File: weather_gridmet_serial.R
# "CHUNKED SERIAL" VERSION
# Optimized for >10,000 points. 
# Processes points in batches to prevent RAM crashes. No Parallel Cluster.

library(terra)
library(sf)
library(dplyr)
library(ncdf4)
library(httr)
library(DSSAT) 
library(lubridate)
library(pbapply) 

process_weather_gridmet <- function(shapefile, start_year, end_year, output_dir, 
                                    id_col, lat_col, lon_col, n_cores, log_file, gridmet_cache_dir,
                                    chunk_size = 3000) { # Added chunk_size argument
  
  message(sprintf("--- Starting GridMET Serial Download (Years: %d-%d) ---", start_year, end_year))
  
  # --- 1. Date Logic ---
  start_date <- as.Date(paste0(start_year, "-01-01"))
  latest_safe_date <- Sys.Date() - 2
  requested_end_date <- as.Date(paste0(end_year, "-12-31"))
  
  if (requested_end_date > latest_safe_date) {
    message(sprintf("NOTICE: Adjusting end date to %s (data availability limit).", latest_safe_date))
    end_date <- latest_safe_date
  } else {
    end_date <- requested_end_date
  }
  
  full_date_seq <- seq(start_date, end_date, by = "day")
  years_needed <- unique(format(full_date_seq, "%Y"))
  
  # --- 2. Setup & Download ---
  gridmet_vars <- list(TMIN = "tmmn", TMAX = "tmmx", RAIN = "pr", SRAD = "srad")
  
  if (!dir.exists(gridmet_cache_dir)) dir.create(gridmet_cache_dir, recursive = TRUE)
  options(timeout = 60000)
  
  downloaded_files <- list()
  for (var_name in names(gridmet_vars)) {
    var_abbrev <- gridmet_vars[[var_name]]
    downloaded_files[[var_name]] <- list()
    for (year in years_needed) {
      file_name <- paste0(var_abbrev, "_", year, ".nc")
      dest_file <- file.path(gridmet_cache_dir, file_name)
      url <- paste0("http://www.northwestknowledge.net/metdata/data/", file_name)
      
      if (!file.exists(dest_file)) {
        tryCatch({ GET(url, write_disk(dest_file, overwrite = TRUE)) }, 
                 error = function(e) message("Download error: ", e$message))
      }
      if (file.exists(dest_file)) downloaded_files[[var_name]][[year]] <- dest_file
    }
  }
  
  # --- 3. PRE-CALCULATE CELL INDICES ---
  message("\n--- Calculating Grid Indices ---")
  
  ref_file <- downloaded_files[["TMIN"]][[years_needed[1]]]
  if(is.null(ref_file)) stop("Fatal: Could not find base TMIN file for spatial reference.")
  
  r_ref <- rast(ref_file)
  
  # Ensure points are in GridMET CRS (EPSG:4326)
  pts_sf <- st_transform(shapefile, 4326)
  coords_mat <- st_coordinates(pts_sf)
  
  # Get the cell numbers for every point
  all_cell_ids <- cellFromXY(r_ref, coords_mat)
  
  # Check for points falling outside the US grid
  valid_mask <- !is.na(all_cell_ids)
  if (sum(!valid_mask) > 0) {
    message(sprintf("WARNING: %d points fall outside the GridMET coverage area. They will be skipped.", sum(!valid_mask)))
  }
  
  rm(r_ref, pts_sf)
  gc()
  
  # --- 4. CHUNKED PROCESSING LOOP ---
  total_points <- nrow(shapefile)
  num_chunks <- ceiling(total_points / chunk_size)
  
  message(sprintf("\n--- Processing %d points in %d chunks (Chunk Size: %d) ---", total_points, num_chunks, chunk_size))
  
  for (k in 1:num_chunks) {
    
    # Define Chunk Indices
    start_idx <- (k - 1) * chunk_size + 1
    end_idx <- min(k * chunk_size, total_points)
    chunk_indices <- start_idx:end_idx
    
    message(sprintf("\nProcessing Chunk %d/%d (Points %d to %d)...", k, num_chunks, start_idx, end_idx))
    
    # Identify Cell IDs for this chunk
    chunk_cell_ids <- all_cell_ids[chunk_indices]
    
    # Skip if all invalid
    if (all(is.na(chunk_cell_ids))) {
      message("Skipping chunk (all points invalid/ocean).")
      next
    }
    
    chunk_data <- list()
    
    # --- EXTRACT DATA FOR THIS CHUNK ---
    for (var_name in names(gridmet_vars)) {
      year_matrices <- list()
      file_list <- downloaded_files[[var_name]]
      sorted_years <- sort(names(file_list))
      
      for (yr in sorted_years) {
        fpath <- file_list[[yr]]
        r <- rast(fpath)
        
        # KEY OPTIMIZATION: Extract only for current CHUNK of cells
        vals <- r[chunk_cell_ids] 
        
        # Handle current year truncation if needed
        days_in_this_year <- seq(as.Date(paste0(yr, "-01-01")), as.Date(paste0(yr, "-12-31")), by="day")
        target_dates_yr <- days_in_this_year[days_in_this_year >= start_date & days_in_this_year <= end_date]
        n_expected <- length(target_dates_yr)
        
        if (ncol(vals) > n_expected) {
          vals <- vals[, 1:n_expected, drop = FALSE]
        }
        
        year_matrices[[yr]] <- as.matrix(vals)
        rm(r, vals); gc() 
      }
      
      chunk_data[[var_name]] <- do.call(cbind, year_matrices)
      rm(year_matrices); gc()
    }
    
    # --- SYNC DATES FOR CHUNK ---
    col_counts <- sapply(chunk_data, ncol)
    min_cols <- min(col_counts)
    chunk_date_seq <- full_date_seq[1:min_cols]
    
    for (v in names(chunk_data)) {
      if (ncol(chunk_data[[v]]) > min_cols) {
        chunk_data[[v]] <- chunk_data[[v]][, 1:min_cols]
      }
    }
    
    # --- WRITE WTH FILES (SERIAL) ---
    message("  -> Writing .WTH files...")
    
    # Standard Serial Loop
    pb <- txtProgressBar(min = 0, max = length(chunk_indices), style = 3)
    
    for (i in 1:length(chunk_indices)) {
      
      global_idx <- chunk_indices[i]
      point_id <- shapefile[[id_col]][global_idx]
      lat <- shapefile[[lat_col]][global_idx]
      lon <- shapefile[[lon_col]][global_idx]
      
      out_f <- file.path(output_dir, paste0(point_id, ".WTH"))
      
      # Skip if exists or invalid
      if (file.exists(out_f) || is.na(chunk_cell_ids[i])) {
        setTxtProgressBar(pb, i)
        next
      }
      
      tryCatch({
        # Grab rows for this point from the chunk matrices
        tmin_row <- chunk_data[["TMIN"]][i, ]
        tmax_row <- chunk_data[["TMAX"]][i, ]
        rain_row <- chunk_data[["RAIN"]][i, ]
        srad_row <- chunk_data[["SRAD"]][i, ]
        
        if (any(is.na(tmin_row))) next
        
        wth_data <- data.frame(
          DATE = chunk_date_seq,
          TMIN = tmin_row - 273.15, 
          TMAX = tmax_row - 273.15,
          RAIN = rain_row,
          SRAD = srad_row * 0.0864  
        )
        
        wth_data$TDEW <- wth_data$TMIN - 2.5
        wth_data$RH2M <- pmin(100, pmax(20, 100 - (wth_data$TMAX - wth_data$TMIN) * 2))
        wth_data$WIND <- -99
        
        tav <- DSSAT::calc_TAV(wth_data)
        amp <- DSSAT::calc_AMP(wth_data)
        wth_data$DATE_FMT <- sprintf("%d%03d", lubridate::year(wth_data$DATE), lubridate::yday(wth_data$DATE))
        
        header <- sprintf(
          "$WEATHER DATA: GridMET Data (Point ID: %s)\n@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n GMET  %8.4f %8.4f   -99 %5.1f %5.1f   -99   -99\n@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
          point_id, lat, lon, tav, amp
        )
        lines <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                         wth_data$DATE_FMT, wth_data$SRAD, wth_data$TMAX, wth_data$TMIN, 
                         wth_data$RAIN, wth_data$TDEW, wth_data$RH2M, wth_data$WIND)
        
        writeLines(c(header, gsub("-99.0", "  -99", lines, fixed=TRUE)), out_f)
        
      }, error = function(e) {
        # Log error but continue
        write(paste("Error on point", point_id, ":", e$message), file = log_file, append = TRUE)
      })
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # --- CLEANUP AFTER CHUNK ---
    rm(chunk_data)
    gc() # Force Garbage Collection
  }
  
  message(sprintf("GridMET processing complete. Output: %s", output_dir))
}
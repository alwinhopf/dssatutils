# ==============================================================================
#  SOIL HELPER: SSURGO (USDA WEB SERVICE) - SMART RESUME & MEMORY OPTIMIZED
#  Filename: soil_ssurgo.R
# ==============================================================================

library(soilDB)
library(sf)
library(dplyr)
library(tidyr)
library(pbapply)
library(parallel)
library(readr)
library(tools)

# --- Wrapper for Retries ---
robust_SDA_query <- function(query, max_retries = 3, retry_delay_seconds = 5) {
  for (attempt in 1:max_retries) {
    result <- try(SDA_query(query), silent = TRUE)
    if (!inherits(result, "try-error")) return(result)
    Sys.sleep(retry_delay_seconds)
  }
  return(NULL)
}

robust_SDA_spatialQuery <- function(point_sf, what, max_retries = 3, retry_delay_seconds = 5) {
  for (attempt in 1:max_retries) {
    result <- try(SDA_spatialQuery(point_sf, what = what), silent = TRUE)
    if (!inherits(result, "try-error")) return(result)
    Sys.sleep(retry_delay_seconds)
  }
  return(NULL)
}

# --- Calculation Logic ---
calculate_soil_properties <- function(soil_properties, top_depth, bottom_depth) {
  ungrouped <- soil_properties %>%
    mutate(
      adj_top = pmax(hzdept_r, top_depth),
      adj_bottom = pmin(hzdepb_r, bottom_depth),
      thickness = adj_bottom - adj_top,
      weighted_clay = claytotal_r * thickness * comppct_r,
      weighted_sand = sandtotal_r * thickness * comppct_r,
      weighted_om = om_r * thickness * comppct_r,
      weighted_bd = dbthirdbar_r * thickness * comppct_r,
      depth_range = paste(top_depth, "-", bottom_depth, "cm")
    ) %>%
    filter(thickness > 0)
  
  grouped <- ungrouped %>%
    group_by(mukey, depth_range) %>%
    summarize(
      clay_pct = sum(weighted_clay, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      sand_pct = sum(weighted_sand, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      silt_pct = 100 - clay_pct - sand_pct,
      om_pct = sum(weighted_om, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      bulk_density = sum(weighted_bd, na.rm = TRUE) / sum(thickness * comppct_r, na.rm = TRUE),
      .groups = 'drop'
    )
  grouped
}

# --- DSSAT Formatting ---
format_dssat_soil_single <- function(profile_data, output_dir) {
  soil_id <- as.character(profile_data$ID[1])
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))
  
  # Redundant check for safety, though main loop handles this now
  if (file.exists(filename)) { return() }
  
  cat("*SOILS: USA SSURGO Soil Profiles\n", file = filename)
  cat("! Generated from SSURGO database\n\n", file = filename, append = TRUE)
  
  cat(sprintf("*%-6s  SSURGO        %9.3f %9.3f\n",
              as.character(soil_id), profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s USA         %9.3f %9.3f \n",
              as.character(soil_id), profile_data$latitude[1], profile_data$longitude[1]),
      file = filename, append = TRUE)
  
  cat("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n", file = filename, append = TRUE)
  cat("    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n", file = filename, append = TRUE)
  
  cat("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n",
      file = filename, append = TRUE)
  
  profile_data %>%
    mutate(depth_num = as.numeric(sub(".*-", "", sub("cm", "", depth_range)))) %>%
    arrange(depth_num) %>%
    group_walk(function(layer, key) {
      depth_val <- layer$depth_num
      slll <- sprintf("%5.3f", layer$SLLL)
      sdul <- sprintf("%5.3f", layer$SDUL)
      ssat <- sprintf("%5.3f", layer$SSAT)
      slll <- sub("^0", " ", slll)
      sdul <- sub("^0", " ", sdul)
      ssat <- sub("^0", " ", ssat)
      
      depth_format <- ifelse(depth_val < 10, sprintf("%6d", depth_val), sprintf("%5d", depth_val))
      
      cat(sprintf("%s   -99 %s %s %s  1.00   -99 %5.2f %5.2f %5.1f %5.1f   -99   -99   -99   -99   -99   -99\n",
                  depth_format, slll, sdul, ssat,
                  layer$bulk_density, layer$om_pct/1.724, layer$clay_pct, layer$silt_pct),
          file = filename, append = TRUE)
    })
  cat("\n", file = filename, append = TRUE)
}

# --- Main Processing Function for SSURGO (Optimized) ---
process_soils_ssurgo <- function(grid_points, output_dir_csv, output_dir_individual, n_cores,
                                 id_col, lat_col, long_col, format_sql_func) {
  
  message("Starting SSURGO Processing (Smart Resume Mode)...")
  
  # --- 1. FILTER: Identify Missing Points BEFORE Loop ---
  # Convert sf to standard dataframe to avoid geometry overhead during filtering
  grid_df <- grid_points %>% 
    st_drop_geometry() %>% 
    as.data.frame()
  
  # Ensure coords are present in dataframe
  if(!lat_col %in% names(grid_df)) grid_df[[lat_col]] <- st_coordinates(grid_points)[,2]
  if(!long_col %in% names(grid_df)) grid_df[[long_col]] <- st_coordinates(grid_points)[,1]
  
  # Get list of existing SOL files (Fast Directory Scan)
  if(dir.exists(output_dir_individual)) {
    existing_files <- tools::file_path_sans_ext(list.files(output_dir_individual, pattern = "\\.SOL$"))
    
    # Identify which points are NOT in existing files
    # Force character conversion to ensure accurate matching
    missing_mask <- ! (as.character(grid_df[[id_col]]) %in% existing_files)
    
    points_to_process <- grid_df[missing_mask, ]
    
    n_total <- nrow(grid_df)
    n_process <- nrow(points_to_process)
    n_skip <- n_total - n_process
    
    message(sprintf("Resume Check: Found %d existing profiles. Processing %d remaining points.", n_skip, n_process))
    
    if (n_process == 0) {
      message("All soil profiles already exist. Skipping SSURGO processing.")
      return(TRUE)
    }
    
  } else {
    # First run, process everything
    points_to_process <- grid_df
    message(sprintf("First Run: Processing all %d points.", nrow(points_to_process)))
  }
  
  # --- Worker Function ---
  process_point_wrapper <- function(point_data_row) {
    ID <- as.character(point_data_row[[id_col]])
    
    # Reconstruct sf object for spatial query
    # We assume WGS84 (EPSG:4326) because previous steps transformed it
    point_sf <- st_as_sf(point_data_row, coords = c(long_col, lat_col), crs = 4326) 
    
    # 1. Spatial Query
    soil_data_query <- robust_SDA_spatialQuery(point_sf, what = 'mukey')
    if (is.null(soil_data_query)) return(NULL)
    
    # 2. Bedrock & Properties
    if (length(soil_data_query$mukey) > 0 && !all(is.na(soil_data_query$mukey))) {
      q_bedrock <- sprintf("SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN %s",
                           format_sql_func(soil_data_query$mukey))
      bedrock_data <- robust_SDA_query(q_bedrock)
      
      bedrock_depth <- 200 # Default
      if (!is.null(bedrock_data)) {
        bd <- as.data.frame(bedrock_data)
        if(nrow(bd) > 0 && !all(is.na(bd$brockdepmin))) bedrock_depth <- min(bd$brockdepmin, na.rm=TRUE)
      }
      if (is.infinite(bedrock_depth)) bedrock_depth <- 200
      
      # Define Layers
      all_layers <- list("0-5cm"=c(0,5), "5-20cm"=c(5,20), "20-35cm"=c(20,35),
                         "35-50cm"=c(35,50), "50-65cm"=c(50,65), "65-80cm"=c(65,80),
                         "80-95cm"=c(80,95), "95-110cm"=c(95,110), "110-125cm"=c(110,125),
                         "125-140cm"=c(125,140), "140-155cm"=c(140,155), "155-170cm"=c(155,170),
                         "170-185cm"=c(170,185), "185-200cm"=c(185,200))
      
      valid_layers <- all_layers[sapply(all_layers, function(x) x[1] < bedrock_depth)]
      if(length(valid_layers) > 0) {
        last <- length(valid_layers)
        if(valid_layers[[last]][2] > bedrock_depth) valid_layers[[last]][2] <- bedrock_depth
        names(valid_layers)[last] <- paste0(valid_layers[[last]][1], "-", valid_layers[[last]][2], "cm")
      } else {
        valid_layers <- list()
        valid_layers[[paste0("0-", bedrock_depth, "cm")]] <- c(0, bedrock_depth)
      }
      
      # Query Properties
      q_soil <- sprintf("SELECT component.mukey, component.cokey, component.comppct_r,
                         chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r,
                         chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r
                         FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey
                         WHERE component.mukey IN %s", format_sql_func(soil_data_query$mukey))
      
      props <- robust_SDA_query(q_soil)
      if(is.null(props)) return(NULL)
      props <- as.data.frame(props)
      if(nrow(props) == 0) return(NULL)
      
      # Calc Props per Layer
      results_list <- lapply(names(valid_layers), function(layer_name) {
        d <- valid_layers[[layer_name]]
        cp <- calculate_soil_properties(props, d[1], d[2])
        if(nrow(cp) > 0) cp$depth_range <- layer_name
        return(cp)
      })
      results_list <- results_list[sapply(results_list, function(x) !is.null(x) && nrow(x) > 0)]
      
      if(length(results_list) > 0) {
        results_df <- do.call(rbind, results_list)
        # Add derived DSSAT physics
        results_df <- results_df %>% mutate(
          ID=ID, longitude=point_data_row[[long_col]], latitude=point_data_row[[lat_col]],
          bedrock_depth_cm=bedrock_depth,
          sand_dec=sand_pct/100, clay_dec=clay_pct/100, om_dec=om_pct/100,
          theta_1500t = -0.024*sand_dec + 0.487*clay_dec + 0.006*om_dec + 0.005*(sand_dec*om_dec) - 0.013*(clay_dec*om_dec) + 0.068*(sand_dec*clay_dec) + 0.031,
          SLLL = theta_1500t + (0.14*theta_1500t - 0.02),
          theta_33t = -0.251*sand_dec + 0.195*clay_dec + 0.011*om_dec + 0.006*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) + 0.452*(sand_dec*clay_dec) + 0.299,
          SDUL = theta_33t + (1.283*(theta_33t)^2 - 0.374*theta_33t - 0.015),
          theta_s33t = 0.278*sand_dec + 0.034*clay_dec + 0.022*om_dec - 0.018*(sand_dec*om_dec) - 0.027*(clay_dec*om_dec) - 0.584*(sand_dec*clay_dec) + 0.078,
          theta_s33 = theta_s33t + (0.636*theta_s33t - 0.107),
          SSAT = SDUL + theta_s33 - 0.097*sand_dec + 0.043
        )
        
        # WRITE .SOL FILE IMMEDIATELY
        format_dssat_soil_single(results_df, output_dir_individual)
        return(results_df)
      }
    }
    return(NULL)
  }
  
  # --- CHUNK SETUP ---
  total_points <- nrow(points_to_process)
  CHUNK_SIZE <- 10000 # Smaller chunk size for better resume granularity
  num_chunks <- ceiling(total_points / CHUNK_SIZE)
  
  message(sprintf("Processing batch of %d points in %d chunks...", total_points, num_chunks))
  
  # --- 1. SETUP CLUSTER (Done ONCE, outside the loop) ---
  if (.Platform$OS.type == "windows") {
    message(sprintf("Initializing cluster with %d cores...", n_cores))
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, { library(soilDB); library(sf); library(dplyr); library(tidyr) })
    clusterExport(cl, varlist=c("process_point_wrapper", "robust_SDA_query", 
                                "robust_SDA_spatialQuery", "calculate_soil_properties",
                                "format_dssat_soil_single", "output_dir_individual", 
                                "id_col", "lat_col", "long_col", "format_sql_func"), 
                  envir = environment())
  } else {
    cl <- n_cores 
  }
  
  # --- 2. CHUNKED LOOP ---
  for (i in 1:num_chunks) {
    start_idx <- (i - 1) * CHUNK_SIZE + 1
    end_idx <- min(i * CHUNK_SIZE, total_points)
    
    message(sprintf("  > Chunk %d/%d (Points %d - %d)", i, num_chunks, start_idx, end_idx))
    
    chunk_data <- points_to_process[start_idx:end_idx, ]
    chunk_list <- split(chunk_data, seq(nrow(chunk_data)))
    
    # Run Parallel on Chunk
    chunk_results <- pblapply(chunk_list, process_point_wrapper, cl = cl)
    
    # Filter NULLs and Combine
    valid_results <- chunk_results[!sapply(chunk_results, is.null)]
    
    if(length(valid_results) > 0) {
      chunk_df <- bind_rows(valid_results)
      
      # Write to CSV (Append Mode)
      write_csv(chunk_df, output_dir_csv, append = file.exists(output_dir_csv))
    }
    
    # --- CRITICAL: CLEAR MEMORY ---
    rm(chunk_list, chunk_results, valid_results)
    if(exists("chunk_df")) rm(chunk_df)
    gc() 
  }
  
  # --- 3. CLEANUP ---
  if (.Platform$OS.type == "windows") {
    stopCluster(cl)
  }
  
  message("SSURGO Processing Complete.")
  return(TRUE) 
}
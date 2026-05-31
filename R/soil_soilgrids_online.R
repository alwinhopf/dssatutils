# ==============================================================================
#  SOIL HELPER: SOILGRIDS ONLINE (REST API & VRT SUPPORT)
#  Filename: soil_soilgrids_online.R
#  Description: Fetches SoilGrids 2.0 data and formats it for DSSAT (.SOL)
# ==============================================================================

# ==============================================================================
#  0. MASTER CONFIGURATION
# ==============================================================================
# USE_REST_API controls the data retrieval strategy inside process_soils_soilgrids_online():
#
#   TRUE  — JSON REST API (https://rest.isric.org)
#           Best for interactive/local use; one request per point; rate-limited.
#           Includes exponential back-off retry logic (up to 5 attempts).
#
#   FALSE — VRT/GDAL virtual rasters (streamed from ISRIC cloud storage)
#           Best for HPC / large batch jobs; avoids per-point HTTP overhead.
#           Requires stable internet; benefits from GDAL block-level caching.
#
# This variable is consumed inside process_soils_soilgrids_online() via:
#   use_rest <- if (exists("USE_REST_API")) USE_REST_API else TRUE
#
# Override from dssat_main_pipeline.R before calling the function:
#   USE_REST_API <- TRUE    # REST API mode  (interactive / local)
#   USE_REST_API <- FALSE   # VRT/GDAL mode  (HPC / batch)
USE_REST_API <- TRUE

# ==============================================================================
#  LIBRARIES
# ==============================================================================
library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(httr)
library(jsonlite)

# ==============================================================================
#  1. PHYSICS & CALCULATIONS (Saxton & Rawls, 2006)
# ==============================================================================
calculate_soil_physics <- function(sand_pct, clay_pct, om_pct) {
  
  # Convert to decimal fractions
  S <- sand_pct / 100
  C <- clay_pct / 100
  OM <- om_pct / 100
  
  # 1. Permanent Wilting Point (LL) / Theta_1500t
  theta_1500t <- -0.024 * S + 0.487 * C + 0.006 * OM + 
    0.005 * (S * OM) - 0.013 * (C * OM) + 
    0.068 * (S * C) + 0.031
  
  # 2. Field Capacity (DUL) / Theta_33t
  theta_33t <- -0.251 * S + 0.195 * C + 0.011 * OM + 
    0.006 * (S * OM) - 0.027 * (C * OM) + 
    0.452 * (S * C) + 0.299
  
  # 3. Saturation (SAT) / Theta_s33t (approximate porosity)
  theta_s33t <- 0.278 * S + 0.034 * C + 0.022 * OM - 
    0.018 * (S * OM) - 0.027 * (C * OM) - 
    0.584 * (S * C) + 0.078
  
  # Adjustments (Rawls canon)
  SLLL <- theta_1500t + (0.14 * theta_1500t - 0.02)
  SDUL <- theta_33t + (1.283 * theta_33t^2 - 0.374 * theta_33t - 0.015)
  #SSAT <- SDUL + theta_s33t - 0.097 * S + 0.043
  #replaced with (Saxton & Rawls 2006) adjustment
  theta_s33 <- theta_s33t + (0.636 * theta_s33t - 0.107)
  SSAT <- SDUL + theta_s33 - 0.097 * S + 0.043
  
  return(list(SLLL = SLLL, SDUL = SDUL, SSAT = SSAT))
}

# ==============================================================================
#  2. DSSAT .SOL FORMATTER (STRICT VALIDATION MODE)
# ==============================================================================
format_dssat_sol_file <- function(site_data, output_dir,
                                  source_name = "ISRIC SoilGrids 2.0",
                                  source_tag = NULL) {

  # 1. Validation: Ensure we actually have data
  if (nrow(site_data) == 0) stop("No soil layers found for this ID.")

  # 2. Strict Check for Critical Missing Values (Sand, Clay, BD, OrgC)
  # DSSAT cannot run if these are NA. We abort this specific file if found.
  critical_vars <- c("sand", "clay", "bdod", "soc_pct")
  if (any(is.na(site_data[, critical_vars]))) {
    stop("Critical soil physical data (Sand, Clay, Bulk Density, or OC) contains NAs.")
  }

  soil_id <- as.character(site_data$ID[1])
  if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  filename <- file.path(output_dir, paste0(soil_id, ".SOL"))

  # --- Header Writing ---
  if (is.null(source_tag))
    source_tag <- ifelse("VRT" %in% names(site_data), "VRT", "REST API")
  cat(sprintf("*SOILS: %s\n", source_name), file = filename)
  cat(sprintf("! Source: %s (%s)\n\n", source_name, source_tag),
      file = filename, append = TRUE)
  
  cat(sprintf("*%-10s  SOILGRIDS     %9.3f %9.3f\n",
              substr(soil_id, 1, 10), site_data$latitude[1], site_data$longitude[1]),
      file = filename, append = TRUE)
  
  cat("@SITE        COUNTRY          LAT     LONG SCS FAMILY\n", file = filename, append = TRUE)
  cat(sprintf(" %-11s World         %9.3f %9.3f \n",
              substr(soil_id, 1, 11), site_data$latitude[1], site_data$longitude[1]),
      file = filename, append = TRUE)
  
  cat("@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE\n", file = filename, append = TRUE)
  cat("    BN   .13     6    .6    73     1     1 IB001 IB001 IB001\n", file = filename, append = TRUE)
  
  cat("@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC\n",
      file = filename, append = TRUE)
  
  # --- Layer Iteration ---
  site_data <- site_data %>% arrange(depth_bottom)
  
  for(i in 1:nrow(site_data)) {
    layer <- site_data[i, ]
    
    # Root Growth Factor
    srgf <- 1.0 * exp(-0.02 * layer$depth_center)
    if (srgf < 0.02) srgf <- 0
    
    # Hydraulic Conductivity (Calculated)
    # Because we validated NAs above, this calculation is now safe.
    ssks <- 60.96 * (10^(0.0126*layer$sand - 0.0064*layer$clay - 0.6))
    if(ssks > 999) ssks <- 999
    
    # Physics check (Ensure physics calc didn't produce NA)
    if(any(is.na(c(layer$SLLL, layer$SDUL, layer$SSAT)))) {
      stop(sprintf("Derived Soil Physics (LL, DUL, SAT) failed for layer %d", i))
    }
    
    cat(sprintf("%6d   -99 %5.3f %5.3f %5.3f %5.2f %5.1f %5.2f %5.2f %5.1f %5.1f %5.1f   -99   -99   -99   -99   -99\n",
                as.integer(layer$depth_bottom),
                layer$SLLL,
                layer$SDUL,
                layer$SSAT,
                srgf,
                ssks,
                layer$bdod,
                layer$soc_pct,
                layer$clay,
                layer$silt,
                layer$cfvo
    ), file = filename, append = TRUE)
  }
  cat("\n", file = filename, append = TRUE)
}

# ==============================================================================
#  3. INTERNAL: REST API WORKER (WITH ROBUST RETRY LOGIC)
# ==============================================================================
fetch_soilgrids_rest <- function(lat, lon) {
  
  url <- "https://rest.isric.org/soilgrids/v2.0/properties/query"
  
  # Flatten parameter list for httr
  base_params <- list(lat = lat, lon = lon, value = "mean")
  props <- c("clay", "sand", "silt", "soc", "bdod", "cfvo")
  depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")
  
  query_params <- base_params
  for(p in props) query_params <- append(query_params, list(property = p))
  for(d in depths) query_params <- append(query_params, list(depth = d))
  
  # --- RETRY CONFIGURATION ---
  max_retries <- 5        # How many times to retry
  base_wait_time <- 2     # Initial wait in seconds
  
  result <- NULL
  
  for(attempt in 1:max_retries) {
    tryCatch({
      
      # Standard polite delay before request
      Sys.sleep(1) 
      
      response <- GET(url, query = query_params)
      status <- status_code(response)
      
      # SUCCESS CASE
      if (status == 200) {
        content_json <- fromJSON(content(response, "text", encoding = "UTF-8"))
        layers_data <- content_json$properties$layers
        
        if(!is.null(layers_data) && length(layers_data) > 0) {
          
          parsed_list <- list()
          
          for(i in 1:nrow(layers_data)) {
            prop_name <- layers_data$name[i]
            depths_df <- layers_data$depths[[i]]
            
            for(j in 1:nrow(depths_df)) {
              depth_label <- depths_df$label[j]
              val <- depths_df$values$mean[j]
              ranges <- as.numeric(unlist(str_extract_all(depth_label, "\\d+")))
              
              parsed_list[[length(parsed_list)+1]] <- data.frame(
                prop = prop_name,
                depth_label = depth_label,
                depth_bottom = ranges[2],
                depth_center = (ranges[1] + ranges[2]) / 2,
                value = val
              )
            }
          }
          result <- bind_rows(parsed_list)
          break # Exit retry loop on success
        }
        
        # RATE LIMIT CASE (429) or SERVER ERROR (503/504)
      } else if (status %in% c(429, 503, 504)) {
        wait_time <- base_wait_time * (2 ^ (attempt - 1)) # Exponential: 2, 4, 8, 16s...
        message(sprintf("Rate limit hit (Status %d). Retrying in %d seconds...", status, wait_time))
        Sys.sleep(wait_time)
        
        # FATAL ERROR CASE (404, 400, etc)
      } else {
        warning(paste("REST API Fatal Error:", status))
        break # Do not retry fatal errors
      }
      
    }, error = function(e) {
      warning(paste("REST Fetch Failed (Network/R Error):", e$message))
      # You might choose to retry here or break depending on preference
      Sys.sleep(2)
    })
  }
  
  return(result)
}

# ==============================================================================
#  4. INTERNAL: VRT WORKER
# ==============================================================================
fetch_soilgrids_vrt <- function(gridfile, id_col) {
  # VRT Setup
  depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm", "60-100cm", "100-200cm")
  depth_centers <- c(2.5, 10, 22.5, 45, 80, 150)
  depth_bottoms <- c(5, 15, 30, 60, 100, 200)
  props <- c("clay", "sand", "silt", "soc", "bdod", "cfvo")
  sg_url <- "/vsicurl/https://files.isric.org/soilgrids/latest/data/"
  
  # Project to IGH
  igh_proj <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
  points_igh <- st_transform(gridfile, igh_proj)
  vect_points <- vect(points_igh)
  
  all_data <- list()
  
  for (prop in props) {
    message(sprintf("VRT: Extracting %s...", prop))
    for (i in seq_along(depths)) {
      depth_label <- depths[i]
      vrt_path <- paste0(sg_url, prop, "/", prop, "_", depth_label, "_mean.vrt")
      
      tryCatch({
        r <- rast(vrt_path)
        vals <- terra::extract(r, vect_points, ID=FALSE)
        
        temp_df <- data.frame(
          ID = gridfile[[id_col]],
          prop = prop,
          depth_label = depth_label,
          depth_bottom = depth_bottoms[i],
          depth_center = depth_centers[i],
          value = vals[,1] 
        )
        all_data[[paste(prop, i, sep="_")]] <- temp_df
      }, error = function(e) message(paste("Skip", prop, depth_label)))
    }
  }
  return(bind_rows(all_data))
}

# ==============================================================================
#  5. MAIN PIPELINE FUNCTION (WITH ERROR LOGGING)
# ==============================================================================
process_soils_soilgrids_online <- function(gridfile, soilfile_csv_path, output_sol_dir, id_col) {
  
  use_rest <- if(exists("USE_REST_API")) USE_REST_API else TRUE
  message(sprintf("--- Starting SoilGrids Extraction (Mode: %s) ---", ifelse(use_rest, "REST API", "VRT")))
  
  if (!inherits(gridfile, "sf")) stop("gridfile must be an sf object")
  
  # WGS84 Coords
  grid_wgs84 <- st_transform(gridfile, 4326)
  coords <- st_coordinates(grid_wgs84)
  grid_wgs84$lon_wgs84 <- coords[,1]
  grid_wgs84$lat_wgs84 <- coords[,2]
  
  full_df <- NULL
  
  # --- BRANCH: REST API (Fetch Loop) ---
  if (use_rest) {
    results_list <- list()
    n_points <- nrow(grid_wgs84)
    
    for(i in 1:n_points) {
      if(i %% 10 == 0) message(sprintf("Processed %d / %d points...", i, n_points))
      
      pid <- grid_wgs84[[id_col]][i]
      # Using the robust fetch function (ensure you updated Part 1/3 previously!)
      res <- fetch_soilgrids_rest(grid_wgs84$lat_wgs84[i], grid_wgs84$lon_wgs84[i])
      
      if(!is.null(res)) {
        res$ID <- pid
        results_list[[i]] <- res
      }
    }
    full_df <- bind_rows(results_list)
    
  } else {
    full_df <- fetch_soilgrids_vrt(gridfile, id_col)
  }
  
  if(is.null(full_df) || nrow(full_df) == 0) {
    stop("No soil data extracted! Check internet connection or coordinates.")
  }
  
  # --- PROCESSING ---
  message("Restructuring and calculating soil physics...")
  
  wide_df <- full_df %>%
    pivot_wider(names_from = prop, values_from = value)
  
  # Conversion & Physics
  processed_df <- wide_df %>%
    mutate(
      clay = clay / 10,       
      sand = sand / 10,       
      silt = silt / 10,       
      soc_pct = soc / 100,    
      om_pct = soc_pct * 1.724,
      bdod = bdod / 100,      
      cfvo = cfvo / 10        
    ) %>%
    rowwise() %>%
    mutate(
      physics = list(calculate_soil_physics(sand, clay, om_pct))
    ) %>%
    unnest_wider(physics) %>%
    ungroup()
  
  coords_df <- data.frame(
    ID = grid_wgs84[[id_col]],
    longitude = grid_wgs84$lon_wgs84,
    latitude = grid_wgs84$lat_wgs84
  )
  
  # Join Lat/Lon back
  # Note: Use inner_join so we only process points that actually retrieved data
  final_df <- inner_join(processed_df, coords_df, by = "ID")
  
  # --- WRITE OUTPUTS (SAFE LOOP) ---
  
  # 1. Master CSV
  mapping_df <- data.frame(ID = grid_wgs84[[id_col]], SOIL_ID = grid_wgs84[[id_col]])
  colnames(mapping_df)[1] <- id_col
  write.csv(mapping_df, soilfile_csv_path, row.names = FALSE)
  
  # 2. Individual SOL files with Error Logging
  unique_ids <- unique(final_df$ID)
  message(sprintf("Writing %d potential .SOL files to: %s", length(unique_ids), output_sol_dir))
  
  # Initialize Log File, save in directory of individual .SOL files
  #log_file_path <- file.path(dirname(output_sol_dir), "soil_processing_errors.log")
  log_file_path <- file.path(output_sol_dir, "soil_processing_errors.log")
  cat(sprintf("Log started: %s\n", Sys.time()), file = log_file_path)
  
  success_count <- 0
  error_count <- 0
  
  for(uid in unique_ids) {
    site_sub <- final_df %>% filter(ID == uid)
    
    # THE SAFE BLOCK
    tryCatch({
      format_dssat_sol_file(site_sub, output_sol_dir)
      success_count <- success_count + 1
    }, error = function(e) {
      # Log the error, but DO NOT STOP the loop
      error_msg <- sprintf("ID: %s | Error: %s", uid, e$message)
      cat(paste0(error_msg, "\n"), file = log_file_path, append = TRUE)
      message(paste("SKIPPED:", error_msg))
      error_count <<- error_count + 1
    })
  }
  
  message(sprintf("SoilGrids processing complete. Success: %d, Errors: %d", success_count, error_count))
  message(sprintf("Check %s for details on skipped points.", log_file_path))
}
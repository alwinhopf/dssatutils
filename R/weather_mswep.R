process_weather_mswep <- function(shapefile, start_year, end_year, output_dir,
                                  id_col, lat_col, lon_col, n_cores, log_file,
                                  mswep_nc_dir) {
  message(sprintf("--- Starting NASA-POWER + MSWEP Processing (Years: %d-%d) ---", start_year, end_year))
  process_weather_nasapower(shapefile, start_year, end_year, output_dir,
                            id_col, lat_col, lon_col, n_cores, log_file)
  pts <- sf::st_transform(shapefile, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  p <- weather_find_nc_file(mswep_nc_dir, c("mswep", "precip", "pr"))
  if (is.na(p)) stop(sprintf("MSWEP precipitation NetCDF not found in %s", mswep_nc_dir))
  rain <- weather_extract_netcdf_series(p, ids, terra::vect(pts), start_year, end_year, "rain")
  for (pid in ids) {
    rf <- rain[[pid]]
    if (is.null(rf) || !length(rf)) next
    f <- file.path(output_dir, paste0(pid, ".WTH"))
    if (!file.exists(f)) next
    ln <- readLines(f, warn = FALSE)
    data_idx <- grep("^\\s*[0-9]{7}\\s+", ln)
    for (ii in data_idx) {
      parts <- strsplit(trimws(ln[ii]), "\\s+")[[1]]
      if (length(parts) >= 5 && parts[1] %in% names(rf)) {
        parts[5] <- sprintf("%.1f", as.numeric(rf[[parts[1]]]))
        ln[ii] <- sprintf("%7s%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f%6.1f",
                          parts[1], as.numeric(parts[2]), as.numeric(parts[3]),
                          as.numeric(parts[4]), as.numeric(parts[5]),
                          as.numeric(parts[6]), as.numeric(parts[7]),
                          as.numeric(parts[8]))
      }
    }
    writeLines(ln, f)
  }
  message(sprintf("\nNASA-POWER + MSWEP processing complete. Check '%s'.\n", output_dir))
}

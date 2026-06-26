.hihydro_first <- function(root, token_groups, depth) {
  for (tokens in token_groups) {
    p <- soil_find_raster(root, tokens, depth)
    if (!is.na(p)) return(p)
  }
  NA_character_
}

process_soils_hihydrosoil <- function(grid_points, hihydrosoil_raster_dir,
                                      output_csv_path, output_sol_dir,
                                      id_col = "ID", lat_col = "LAT",
                                      long_col = "LONG", integer_scale = 0.0001) {
  if (!dir.exists(hihydrosoil_raster_dir))
    stop(sprintf("HiHydroSoil raster directory not found: %s", hihydrosoil_raster_dir))
  pts <- sf::st_transform(grid_points, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)
  rows <- list()
  for (j in seq_len(nrow(SOIL_RASTER_DEPTHS))) {
    d <- SOIL_RASTER_DEPTHS[j, ]
    slll <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("pf42", "pf4.2", "pF4.2"), d$token), pts_vect, integer_scale)
    sdul <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("pf25", "pf2.5", "pF2.5", "pf2", "pF2"), d$token), pts_vect, integer_scale)
    ssat <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("thetas", "theta_s", "sat"), d$token), pts_vect, integer_scale)
    ksat <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("ksat", "conductivity"), d$token), pts_vect, integer_scale)
    om <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("organic", "om"), d$token), pts_vect, integer_scale)
    tex <- soil_sample_raster(.hihydro_first(hihydrosoil_raster_dir, list("texture", "usda"), d$token), pts_vect)
    sand <- soil_sample_raster(soil_find_raster(hihydrosoil_raster_dir, "sand", d$token), pts_vect)
    clay <- soil_sample_raster(soil_find_raster(hihydrosoil_raster_dir, "clay", d$token), pts_vect)
    silt <- soil_sample_raster(soil_find_raster(hihydrosoil_raster_dir, "silt", d$token), pts_vect)
    for (i in seq_along(ids)) {
      sa <- sand[i]; si <- silt[i]; cl <- clay[i]
      if (!is.finite(sa) || !is.finite(si) || !is.finite(cl)) {
        pct <- soil_texture_to_pct(tex[i]); sa <- pct[1]; si <- pct[2]; cl <- pct[3]
      }
      rows[[length(rows) + 1]] <- data.frame(
        ID = ids[i], latitude = lats[i], longitude = lons[i],
        depth_bottom = d$bottom, depth_center = d$center,
        sand = sa, clay = cl, silt = si,
        bdod = ifelse(is.finite(ssat[i]), (1 - ssat[i]) * 2.65, NA_real_),
        soc_pct = ifelse(is.finite(om[i]), om[i] / 1.724, 1.0),
        cfvo = 0, SLLL = slll[i], SDUL = sdul[i], SSAT = ssat[i],
        # ksat already scaled to cm/day; DSSAT SSKS is cm/h (/24).
        SSKS = ifelse(is.finite(ksat[i]), ksat[i] / 24.0, NA_real_))
    }
  }
  df <- do.call(rbind, rows)
  df <- df[complete.cases(df[, c("sand", "clay", "bdod", "SLLL", "SDUL", "SSAT")]), ]
  if (!nrow(df)) stop("No usable HiHydroSoil data extracted. Check raster names, scaling, and coordinates.")
  df$SDUL <- pmax(df$SDUL, df$SLLL + 0.04)
  df$SSAT <- pmax(df$SSAT, df$SDUL + 0.04)
  soil_write_mapping(ids, output_csv_path)
  soil_write_profiles(df, output_sol_dir, "HiHydroSoil v2.0", "local hydraulic rasters")
}

# Internal helpers for local-raster soil sources.

SOIL_RASTER_DEPTHS <- data.frame(
  token = c("0_5", "5_15", "15_30", "30_60", "60_100", "100_200"),
  bottom = c(5, 15, 30, 60, 100, 200),
  center = c(2.5, 10, 22.5, 45, 80, 150)
)

soil_find_raster <- function(root, tokens, depth_token = NULL) {
  files <- list.files(root, pattern = "\\.(tif|tiff|vrt)$", full.names = TRUE,
                      recursive = TRUE, ignore.case = TRUE)
  base <- tolower(basename(files))
  hit <- Reduce(`&`, lapply(tolower(tokens), function(t) grepl(t, base, fixed = TRUE)))
  if (!is.null(depth_token)) {
    d1 <- tolower(depth_token); d2 <- gsub("_", "-", d1, fixed = TRUE)
    hit <- hit & (grepl(d1, base, fixed = TRUE) | grepl(d2, base, fixed = TRUE))
  }
  if (any(hit)) files[which(hit)[1]] else NA_character_
}

soil_sample_raster <- function(path, pts_vect, scale = 1) {
  if (is.na(path) || !file.exists(path)) return(rep(NA_real_, nrow(pts_vect)))
  vals <- tryCatch(as.numeric(terra::extract(terra::rast(path), pts_vect, ID = FALSE)[, 1]),
                   error = function(e) rep(NA_real_, nrow(pts_vect)))
  vals * scale
}

soil_texture_to_pct <- function(cls) {
  table <- rbind(c(92,5,3), c(82,12,6), c(65,25,10), c(43,39,18),
                 c(20,65,15), c(10,80,10), c(52,7,41), c(45,15,40),
                 c(32,34,34), c(20,40,40), c(10,34,56), c(22,20,58))
  if (is.na(cls)) return(c(NA_real_, NA_real_, NA_real_))
  idx <- as.integer(round(cls))
  if (idx < 1 || idx > nrow(table)) return(c(NA_real_, NA_real_, NA_real_))
  table[idx, ]
}

soil_write_mapping <- function(ids, output_csv_path) {
  utils::write.csv(data.frame(ID = ids, SOIL_ID = ids), output_csv_path, row.names = FALSE)
}

soil_write_profiles <- function(df, output_sol_dir, source_name, source_tag) {
  for (uid in unique(df$ID)) {
    format_dssat_sol_file(df[df$ID == uid, , drop = FALSE], output_sol_dir,
                          source_name = source_name, source_tag = source_tag)
  }
}

soil_add_physics <- function(df) {
  phys <- lapply(seq_len(nrow(df)), function(i)
    calculate_soil_physics(df$sand[i], df$clay[i], df$soc_pct[i] * 1.724))
  df$SLLL <- vapply(phys, `[[`, numeric(1), "SLLL")
  df$SDUL <- vapply(phys, `[[`, numeric(1), "SDUL")
  df$SSAT <- vapply(phys, `[[`, numeric(1), "SSAT")
  df
}

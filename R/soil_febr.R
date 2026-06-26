# Soil source: FEBR / Embrapa Brazilian soil profiles (processed CSV).
# Brazil-specific soil pairing with BR-DWGD/Xavier weather; nearest-profile
# ingestion like WoSIS. Required CSV columns: profile_id, latitude, longitude,
# depth_bottom, sand, clay, silt, bdod, soc_pct.

process_soils_febr <- function(grid_points, febr_profile_csv, output_csv_path,
                               output_sol_dir, id_col = "ID", lat_col = "LAT",
                               long_col = "LONG") {
  process_soils_wosis(grid_points, febr_profile_csv, output_csv_path, output_sol_dir,
                      id_col = id_col, lat_col = lat_col, long_col = long_col)
}

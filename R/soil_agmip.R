# ==============================================================================
#  SOIL HELPER: AGMIP / HAN GLOBAL DSSAT .SOL -> per-point .SOL + mapping CSV
# ==============================================================================

#' Process the AgMIP / Han global DSSAT soil profile database
#'
#' The AgMIP/Han source is distributed as country-level DSSAT-ready .SOL files
#' at 5 arc-min (~10 km) resolution. Because the profiles already include DSSAT
#' hydraulic properties, this wrapper delegates to the generic external-.SOL
#' mapper used by process_soils_soilgrids().
#'
#' Data DOI: https://doi.org/10.7910/DVN/1PEEY0
#' Paper: Han et al. 2019, Environ. Model. Softw. 119:70-83
process_soils_agmip <- function(
  grid_points,
  source_sol_file,
  output_csv_path,
  output_sol_dir,
  id_col = "ID",
  numeric_only_ids = TRUE,
  numeric_width = 8
) {
  message("--- Starting AgMIP/Han DSSAT Soil Profile Extraction ---")
  process_soils_soilgrids(
    grid_points = grid_points,
    source_sol_file = source_sol_file,
    output_csv_path = output_csv_path,
    output_sol_dir = output_sol_dir,
    id_col = id_col,
    numeric_only_ids = numeric_only_ids,
    numeric_width = numeric_width
  )
}

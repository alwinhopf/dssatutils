# File: soil_agmip.py
# ---------------------------------------------------------------------------
# Soil source: AgMIP / Han et al. global DSSAT-ready soil profile database.
#
# The source data are distributed as country-level DSSAT .SOL files at 5 arc-min
# (~10 km) resolution. They already contain DSSAT hydraulic properties to 2 m,
# so this backend intentionally reuses the generic external-.SOL nearest-profile
# mapper instead of recomputing soil physics.
#
# Data DOI: https://doi.org/10.7910/DVN/1PEEY0
# Paper:    Han et al. 2019, Environ. Model. Softw. 119:70-83
# ---------------------------------------------------------------------------

from .soil_soilgrids import process_soils_soilgrids


def process_soils_agmip(
    grid_points,
    source_sol_file: str,
    output_csv_path: str,
    output_sol_dir: str,
    id_col: str = "ID",
    numeric_only_ids: bool = True,
    numeric_width: int = 8,
):
    """Build per-point DSSAT .SOL files from the AgMIP/Han global profile DB.

    This is a semantic wrapper around :func:`process_soils_soilgrids`: both data
    products are preformatted DSSAT master ``.SOL`` files, and the correct
    operation for a gridded run is to select the nearest source profile for each
    simulation point, rewrite its soil ID, and emit the standard mapping CSV.
    """
    print("--- Starting AgMIP/Han DSSAT Soil Profile Extraction ---")
    return process_soils_soilgrids(
        grid_points=grid_points,
        source_sol_file=source_sol_file,
        output_csv_path=output_csv_path,
        output_sol_dir=output_sol_dir,
        id_col=id_col,
        numeric_only_ids=numeric_only_ids,
        numeric_width=numeric_width,
    )

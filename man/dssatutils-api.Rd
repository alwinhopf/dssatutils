\name{dssatutils-api}
\alias{dssatutils-api}
\alias{audit_weather_file_quality}
\alias{audit_weather_quality}
\alias{era5land_set_cds_key}
\alias{extract_chirps_v3_rainfall}
\alias{merge_rainfall_into_weather}
\alias{is_wth_valid}
\alias{process_soils_agmip}
\alias{process_soils_china}
\alias{process_soils_esdb}
\alias{process_soils_febr}
\alias{process_soils_gnatsgo}
\alias{process_soils_gsde}
\alias{process_soils_hihydrosoil}
\alias{process_soils_hwsd}
\alias{process_soils_isdasoil}
\alias{process_soils_lucas}
\alias{process_soils_openlandmap}
\alias{process_soils_polaris}
\alias{process_soils_slc}
\alias{process_soils_slga}
\alias{process_soils_soilgrids}
\alias{process_soils_soilgrids_online}
\alias{process_soils_ssurgo}
\alias{process_soils_ssurgo_alderman}
\alias{process_soils_wise30sec}
\alias{process_soils_wosis}
\alias{process_weather_agcfsr}
\alias{process_weather_agera5}
\alias{process_weather_agmerra}
\alias{process_weather_anusplin}
\alias{process_weather_aphrodite}
\alias{process_weather_chelsa_w5e5}
\alias{process_weather_cmfd}
\alias{process_weather_crujra}
\alias{process_weather_daymet}
\alias{process_weather_dwd}
\alias{process_weather_eobs}
\alias{process_weather_era5_land}
\alias{process_weather_ghcn}
\alias{process_weather_gridmet}
\alias{process_weather_merra2}
\alias{process_weather_mswep}
\alias{process_weather_mswx}
\alias{process_weather_nasapower}
\alias{process_weather_nasapower_chirps}
\alias{process_weather_nasapower_chirps_v3}
\alias{process_weather_openmeteo}
\alias{process_weather_pgf}
\alias{process_weather_prism}
\alias{process_weather_silo}
\alias{process_weather_tamsat}
\alias{process_weather_terraclimate}
\alias{process_weather_xavier}
\alias{pull_profile_by_coords_alderman}
\alias{pull_profile_by_name_alderman}
\alias{repair_weather_date_gaps}
\alias{repair_weather_file_date_gaps}
\alias{repair_weather_file_missing_values}
\alias{repair_weather_file_temperature_inversions}
\alias{repair_weather_missing_values}
\alias{repair_weather_temperature_inversions}
\alias{setup_cds_credentials}
\title{DSSAT Weather and Soil Utility API}
\description{
Public weather and soil processors, quality checks, repair helpers, and
credential setup functions. Each processor writes DSSAT-compatible weather or
soil files; source-specific arguments and examples are documented in the
package README.
}
\details{
Weather processors accept point coordinates, a requested calendar period, an
output directory, column names, worker count, and a log path. They write one
complete DSSAT \code{.WTH} file per point and fail when required TMAX, TMIN,
RAIN, or SRAD forcing is incomplete. Soil processors likewise write DSSAT
\code{.SOL} profiles and mapping tables. Network credentials and source-specific
local inputs are never inferred from another provider.

R and Python expose the same public source names. See the package README data
source matrix for coverage, credentials, optional dependencies, cache behavior,
and complete call examples.
}
\value{
Processors write DSSAT input files and return source-specific metadata or an
invisible completion value. Quality-audit functions return structured issue
summaries; credential setup returns the effective CDS credential path and URL.
}
\examples{
\dontrun{
setup_cds_credentials(prompt = TRUE)
process_weather_openmeteo(points, 2001, 2002, "weather",
  "ID", "LAT", "LONG", 2, "weather_errors.log")
}
}

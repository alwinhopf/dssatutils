# Symbols evaluated within dplyr/data.table expressions. Declaring them keeps
# package checks meaningful without changing their run-time lookup semantics.
utils::globalVariables(c(
  ":=", ".", ".data", "adj_bottom", "adj_top", "albedodry_r",
  "ALLSKY_SFC_SW_DWN", "bdod", "bedrock", "bulk_density", "cfvo",
  "CHIRPS_RESOLUTION", "clay", "clay_dec", "clay_pct", "claytotal_r",
  "coarse_fraction", "cokey", "comppct_r", "DATE_obj", "dbfifteenbar_r",
  "dbovendry_r", "dbtenthbar_r", "dbthirdbar_r", "depth_bottom", "depth_num",
  "depth_range", "DEW_C", "DOY", "fragvol_r", "fragvol_raw", "hydgrp",
  "hzdepb_r", "hzdept_r", "hzname", "i", "ID", "ksat_r", "lat",
  "latitude", "layer_tbl", "lon", "longitude", "MM", "mukey", "om_dec",
  "om_pct", "om_r", "partdensity", "physics", "PRECTOTCORR", "prop",
  "RAIN", "RAIN_MM", "RH2M", "SADC", "sand", "sand_dec", "sand_pct",
  "sandtotal_r", "SBDM", "SCEC", "SDUL", "SDUL_ptf", "SDUL_raw", "silt",
  "silt_pct", "silttotal_r", "SLB", "SLCF", "SLCL", "SLHB", "SLHW",
  "SLLL", "SLLL_ptf", "SLLL_raw", "SLMH", "SLNI", "SLOC", "slope_r",
  "SLSI", "soc", "soc_pct", "SOIL_ID", "SOIL_LAT", "SOIL_LON",
  "SOURCE_SOIL_ID", "SRAD", "SRAD_MJ", "SRGF", "SSAT", "SSAT_ptf",
  "SSAT_raw", "SSKS", "T2M_C", "T2M_MAX", "T2M_MIN", "T2MDEW", "TAVG",
  "TAVG_MON", "TDEW", "theta_1500t", "theta_33t", "theta_s33",
  "theta_s33t", "thickness", "TMAX", "TMIN", "value", "weighted_bd",
  "weighted_clay", "weighted_om", "weighted_sand", "wfifteenbar_r", "WIND",
  "WIND_MS", "WS2M", "wsatiated_r", "wtenthbar_r", "wthirdbar_r", "YEAR"
))

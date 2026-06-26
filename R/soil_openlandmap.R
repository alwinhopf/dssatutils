# Soil source: OpenLandMap / OpenGeoHub global soil COGs (live remote sampling).
# Samples cloud-optimized GeoTIFFs over HTTP via the OpenLandMap STAC catalog
# (URLs resolved at runtime). Texture from the 120 m global_soil_props product;
# bulk density + organic carbon from the 250 m COGs at the nearest point-depth.

OLM_STAC_BASE <- "https://s3.eu-central-1.wasabisys.com/stac/openlandmap/"
OLM_COLLECTIONS <- c(
  clay = "clay.tot_iso.11277.2020.wpct",
  sand = "sand.tot_iso.11277.2020.wpct",
  silt = "silt.tot_iso.11277.2020.wpct",
  bdod = "bulkdens.fineearth_usda.4a1h",
  oc   = "organic.carbon_usda.6a1c"
)
OLM_LAYERS <- data.frame(bottom = c(30, 60, 100), center = c(15, 45, 80))

.olm_depth_mid <- function(key) {
  m <- regmatches(key, regexec("b([0-9]+)cm\\.\\.([0-9]+)cm", key))[[1]]
  if (length(m) == 3) return((as.numeric(m[2]) + as.numeric(m[3])) / 2)
  m <- regmatches(key, regexec("b([0-9]+)cm", key))[[1]]
  if (length(m) == 2) return(as.numeric(m[2]))
  NA_real_
}

.olm_resolve_assets <- function(collection) {
  col <- jsonlite::fromJSON(paste0(OLM_STAC_BASE, collection, "/collection.json"), simplifyVector = FALSE)
  item_link <- Filter(function(l) identical(l$rel, "item"), col$links)[[1]]$href
  item_url <- paste0(OLM_STAC_BASE, collection, "/", sub("^\\./", "", item_link))
  assets <- jsonlite::fromJSON(item_url, simplifyVector = FALSE)$assets
  best <- list()
  for (k in names(assets)) {
    href <- assets[[k]]$href
    if (grepl("_m_", k) && grepl("\\.tif$", href) && !grepl("preview", k)) {
      mid <- .olm_depth_mid(k)
      if (!is.na(mid)) {
        fine <- grepl("120m", k)
        key <- as.character(mid)
        if (is.null(best[[key]]) || (fine && !best[[key]]$fine))
          best[[key]] <- list(mid = mid, href = href, fine = fine)
      }
    }
  }
  do.call(rbind, lapply(best, function(b) data.frame(mid = b$mid, href = b$href, stringsAsFactors = FALSE)))
}

.olm_nearest <- function(assets, target) assets$href[which.min(abs(assets$mid - target))]

.olm_sample <- function(url, pts_vect) {
  r <- tryCatch(terra::rast(paste0("/vsicurl/", url)), error = function(e) NULL)
  if (is.null(r)) return(rep(NA_real_, nrow(pts_vect)))
  v <- terra::extract(r, pts_vect, ID = FALSE)[, 1]
  nd <- terra::NAflag(r)
  if (!is.na(nd)) v[v == nd] <- NA_real_
  as.numeric(v)
}

process_soils_openlandmap <- function(grid_points, output_csv_path, output_sol_dir,
                                      id_col = "ID", lat_col = "LAT", long_col = "LONG") {
  message("--- Starting OpenLandMap COG Extraction (live) ---")
  pts <- sf::st_transform(grid_points, 4326)
  ids <- as.character(sf::st_drop_geometry(pts)[[id_col]])
  xy <- sf::st_coordinates(pts); lats <- xy[, 2]; lons <- xy[, 1]
  pts_vect <- terra::vect(pts)
  assets <- lapply(OLM_COLLECTIONS, .olm_resolve_assets)

  rows <- list()
  for (j in seq_len(nrow(OLM_LAYERS))) {
    dctr <- OLM_LAYERS$center[j]; dbot <- OLM_LAYERS$bottom[j]
    clay <- .olm_sample(.olm_nearest(assets$clay, dctr), pts_vect)
    sand <- .olm_sample(.olm_nearest(assets$sand, dctr), pts_vect)
    silt <- .olm_sample(.olm_nearest(assets$silt, dctr), pts_vect)
    bd <- .olm_sample(.olm_nearest(assets$bdod, dctr), pts_vect)
    bd <- ifelse(bd > 10, bd / 100, bd)
    oc <- .olm_sample(.olm_nearest(assets$oc, dctr), pts_vect)
    soc_pct <- oc / 50  # x5 g/kg -> g/kg (/5) -> percent (/10)
    silt <- ifelse(!is.finite(silt) & is.finite(sand) & is.finite(clay),
                   pmax(0, 100 - sand - clay), silt)
    rows[[j]] <- data.frame(ID = ids, latitude = lats, longitude = lons,
                            depth_bottom = dbot, depth_center = dctr,
                            sand = sand, clay = clay, silt = silt,
                            bdod = bd, soc_pct = soc_pct, cfvo = 0)
  }
  df <- do.call(rbind, rows)
  df <- df[complete.cases(df[, c("sand", "clay", "bdod", "soc_pct")]), ]
  if (!nrow(df)) stop("No usable OpenLandMap data sampled. Check connectivity / coordinates.")
  df <- soil_add_physics(df)
  soil_write_mapping(ids, output_csv_path)
  soil_write_profiles(df, output_sol_dir, "OpenLandMap (OpenGeoHub)", "STAC COG sampling")
}

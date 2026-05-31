#!/usr/bin/env bash
# PHASE 2 — rewire the "DSSAT ML Phenology Prediction" repo to consume dssatutils.
# Idempotent and conservative: backs up touched files, removes only the duplicated
# weather/soil scripts, and inserts a guarded install+library() in pipeline/config.R.
#
# Run AFTER the dssatutils remote+tag exist (or set USE_LOCAL=1 to install from path).
set -uo pipefail

REPO="/Users/alwinhopf/Documents/GitHub/DSSAT ML Phenology Prediction"
PKG_TAG="alwinhopf/dssatutils@v0.1.0"
USE_LOCAL="${USE_LOCAL:-0}"          # 1 = install from local path instead of GitHub
LOCAL_PKG="/Users/alwinhopf/Documents/GitHub/dssatutils"
STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$REPO/.migration_backup_$STAMP"

cd "$REPO" || { echo "repo not found: $REPO"; exit 1; }
mkdir -p "$BACKUP"
echo "Backup dir: $BACKUP"

# 1) Remove duplicated weather/soil R scripts (back them up first).
for f in r_scripts/weather_daymet.R r_scripts/weather_gridmet.R \
         r_scripts/weather_nasapower.R r_scripts/soil_soilgrids.R \
         r_scripts/soil_soilgrids_online.R r_scripts/soil_ssurgo.R; do
  if [ -f "$f" ]; then mkdir -p "$BACKUP/$(dirname "$f")"; cp "$f" "$BACKUP/$f"; git rm -q "$f" 2>/dev/null || rm -f "$f"; echo "removed $f"; fi
done
# NOTE: this repo also has weather_agera5/openmeteo/chirps & soil_hwsd ONLY if it
# was synced from the Gridded repo. Remove any that exist:
for f in r_scripts/weather_agera5.R r_scripts/weather_openmeteo.R \
         r_scripts/weather_nasapower_chirps.R r_scripts/soil_hwsd.R; do
  if [ -f "$f" ]; then cp "$f" "$BACKUP/$f"; git rm -q "$f" 2>/dev/null || rm -f "$f"; echo "removed $f"; fi
done

# 2) Insert guarded loader into pipeline/config.R (only once).
CFG="pipeline/config.R"
if [ -f "$CFG" ] && ! grep -q "dssatutils" "$CFG"; then
  cp "$CFG" "$BACKUP/config.R"
  if [ "$USE_LOCAL" = "1" ]; then
    INSTALL_LINE="  remotes::install_local(\"$LOCAL_PKG\", force = FALSE)"
  else
    INSTALL_LINE="  remotes::install_github(\"$PKG_TAG\")"
  fi
  cat > /tmp/dssatutils_loader.R <<EOF
# --- Shared weather/soil utilities (extracted to dssatutils) ---
if (!requireNamespace("dssatutils", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
$INSTALL_LINE
}
library(dssatutils)
# --- end dssatutils loader ---
EOF
  # Prepend loader to top of config.R
  cat /tmp/dssatutils_loader.R "$CFG" > /tmp/config_new.R && mv /tmp/config_new.R "$CFG"
  echo "patched $CFG with dssatutils loader"
else
  echo "config.R missing or already patched — skipping"
fi

# 3) Replace any source("r_scripts/weather_*|soil_*") with nothing (function now from pkg).
#    Conservative: comment them out rather than delete, so context is preserved.
grep -rln 'source(.*r_scripts/\(weather\|soil\)_' --include=*.R . 2>/dev/null | while read -r file; do
  cp "$file" "$BACKUP/$(basename "$file").bak"
  sed -i '' -E 's/^([[:space:]]*)(source\(.*r_scripts\/(weather|soil)_.*)/\1# [dssatutils] \2/' "$file"
  echo "commented source() lines in $file"
done

echo
echo "PHASE 2 edits staged. Review with:  git -C \"$REPO\" diff"
echo "Then in R:  renv::init(); renv::snapshot()   # repo had no lockfile"
echo "Backups in: $BACKUP"

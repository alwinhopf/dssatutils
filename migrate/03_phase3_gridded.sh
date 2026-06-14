#!/usr/bin/env bash
# PHASE 3 — rewire the "DSSAT Gridded Run Tutorial" repo (source of truth) to
# consume dssatutils instead of its local r_scripts/ + python_scripts/ copies.
# Landcover scripts STAY. Idempotent; backs up everything it touches.
#
# Run AFTER dssatutils remote+tag exist (or USE_LOCAL=1 for local path install).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE_ROOT="${WORKSPACE_ROOT:-${GITHUB_WORKSPACE:-$(cd "$SCRIPT_DIR/../.." && pwd)}}"
REPO="${REPO:-$WORKSPACE_ROOT/DSSAT_Gridded_Run_Tutorial}"
PKG_TAG="alwinhopf/dssatutils@v0.1.0"
PIP_PIN="dssatutils @ git+https://github.com/alwinhopf/dssatutils.git@v0.1.0"
USE_LOCAL="${USE_LOCAL:-0}"
LOCAL_PKG="${LOCAL_PKG:-$WORKSPACE_ROOT/dssatutils}"
STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$REPO/.migration_backup_$STAMP"

cd "$REPO" || { echo "repo not found: $REPO"; exit 1; }
mkdir -p "$BACKUP"
echo "Backup dir: $BACKUP"

WEATHER_SOIL=(weather_daymet weather_gridmet weather_nasapower weather_agera5 \
  weather_openmeteo weather_nasapower_chirps soil_soilgrids soil_soilgrids_online \
  soil_ssurgo soil_hwsd)

# ---- R pipeline: replace the 10 weather/soil source() lines with library() ----
RPIPE="dssat_main_pipeline.R"
if [ -f "$RPIPE" ]; then
  cp "$RPIPE" "$BACKUP/$(basename "$RPIPE")"
  # Comment out weather_/soil_ source() lines (KEEP landcover_*).
  sed -i '' -E 's/^([[:space:]]*)(source\(file\.path\(SCRIPT_DIR, "(weather|soil)_[a-z0-9_]*\.R"\)\).*)/\1# [dssatutils] \2/' "$RPIPE"
  # Insert library(dssatutils) once, just before the first commented line.
  if ! grep -q 'library(dssatutils)' "$RPIPE"; then
    awk '!ins && /# \[dssatutils\] source\(file\.path\(SCRIPT_DIR, "(weather|soil)_/ {print "library(dssatutils)  # [dssatutils] shared weather/soil sources"; ins=1} {print}' "$RPIPE" > /tmp/rpipe_new && mv /tmp/rpipe_new "$RPIPE"
  fi
  echo "patched $RPIPE (landcover source() lines left intact)"
fi

# ---- Python pipeline: swap flat imports for package imports --------------------
PYPIPE="dssat_main_pipeline.py"
if [ -f "$PYPIPE" ]; then
  cp "$PYPIPE" "$BACKUP/$(basename "$PYPIPE")"
  for m in "${WEATHER_SOIL[@]}"; do
    # from <m> import process_xxx   ->  from dssatutils.<m> import process_xxx
    sed -i '' -E "s/^from $m([[:space:]]+)import/from dssatutils.$m\1import/" "$PYPIPE"
  done
  # lazy/local imports deeper in the file
  sed -i '' -E 's/^([[:space:]]*)import soil_soilgrids_online as _sg_mod/\1import dssatutils.soil_soilgrids_online as _sg_mod/' "$PYPIPE"
  sed -i '' -E 's/^([[:space:]]*)import weather_nasapower_chirps as _wc/\1import dssatutils.weather_nasapower_chirps as _wc/' "$PYPIPE"
  echo "patched $PYPIPE imports"
fi

# ---- Remove the moved source files (KEEP landcover_*) -------------------------
for m in "${WEATHER_SOIL[@]}"; do
  for ext in R py; do
    sub=$([ "$ext" = R ] && echo r_scripts || echo python_scripts)
    f="$sub/$m.$ext"
    if [ -f "$f" ]; then mkdir -p "$BACKUP/$sub"; cp "$f" "$BACKUP/$f"; git rm -q "$f" 2>/dev/null || rm -f "$f"; echo "removed $f"; fi
  done
done
echo "kept: r_scripts/landcover_*  python_scripts/landcover_*"

# ---- Dependency pins ----------------------------------------------------------
if [ -f requirements.txt ] && ! grep -q 'dssatutils' requirements.txt; then
  cp requirements.txt "$BACKUP/requirements.txt"
  if [ "$USE_LOCAL" = "1" ]; then
    echo "-e $LOCAL_PKG" >> requirements.txt
  else
    echo "$PIP_PIN" >> requirements.txt
  fi
  echo "added dssatutils pin to requirements.txt"
fi
if [ -f setup_renv.R ] && ! grep -q 'dssatutils' setup_renv.R; then
  cp setup_renv.R "$BACKUP/setup_renv.R"
  if [ "$USE_LOCAL" = "1" ]; then
    echo 'remotes::install_local("'"$LOCAL_PKG"'", force = FALSE)' >> setup_renv.R
  else
    echo 'remotes::install_github("'"$PKG_TAG"'")' >> setup_renv.R
  fi
  echo "added dssatutils install to setup_renv.R"
fi

echo
echo "PHASE 3 edits staged. Review with:  git -C \"$REPO\" diff"
echo "Backups in: $BACKUP"

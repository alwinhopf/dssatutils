#!/usr/bin/env bash
# PHASE 2 — rewire "DSSAT ML Phenology Prediction" to consume dssatutils.
#
# This repo does NOT use literal `source("r_scripts/...")`. It sources via a
# variable, in several files, and sometimes with `local = TRUE`:
#   source(file.path(R_SCRIPTS_DIR, "soil_soilgrids_online.R"))          # 01_particle_filter.R, 04_cohesive_calibration.R
#   source(file.path(R_SCRIPTS_DIR, "soil_soilgrids.R"), local = TRUE)   # utils.R
#   source(file.path(R_SCRIPTS_DIR, "soil_soilgrids_online.R"), local = TRUE) # utils.R
#   source(file.path(PROJECT_ROOT, "r_scripts/soil_soilgrids_online.R")) # scratch/run_g29_debug.R
#
# Strategy: replace every weather_/soil_ source() line (any of those forms) with
# `suppressMessages(library(dssatutils))` — preserving indentation — so the
# functions resolve at each call site regardless of the config-sourcing chain
# (library() attaches globally, which is also fine inside a `local = TRUE` fn).
# Then prepend a one-time guarded install+library to pipeline/config.R, and
# delete the now-duplicated weather/soil files from r_scripts/ (landcover stays).
#
# Idempotent; backs up every touched file. Run from anywhere.
set -uo pipefail

REPO="/Users/alwinhopf/Documents/GitHub/DSSAT ML Phenology Prediction"
PKG_TAG="alwinhopf/dssatutils@v0.1.0"
USE_LOCAL="${USE_LOCAL:-0}"                 # 1 = install from local path (no remote yet)
LOCAL_PKG="/Users/alwinhopf/Documents/GitHub/dssatutils"
STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="$REPO/.migration_backup_$STAMP"

cd "$REPO" || { echo "repo not found: $REPO"; exit 1; }
mkdir -p "$BACKUP"
echo "Backup dir: $BACKUP"

# 1) Rewrite weather_/soil_ source() lines -> suppressMessages(library(dssatutils)).
#    Matches both the R_SCRIPTS_DIR form and the literal r_scripts/ form, with or
#    without a trailing `, local = TRUE)`. Keeps leading indentation.
# (while-read loop, not mapfile — macOS /bin/bash is 3.2 and has no mapfile)
grep -rlE 'source\([^)]*(R_SCRIPTS_DIR, "|r_scripts/)(weather|soil)_[a-z0-9_]*\.R' --include="*.R" . 2>/dev/null | while IFS= read -r file; do
  rel="${file#./}"
  mkdir -p "$BACKUP/$(dirname "$rel")"
  cp "$file" "$BACKUP/$rel"
  # Replace any matching whole source(...) line with the library() call.
  sed -i '' -E 's/^([[:space:]]*)source\([^)]*(R_SCRIPTS_DIR, "|r_scripts\/)(weather|soil)_[a-z0-9_]*\.R"?[^)]*\).*/\1suppressMessages(library(dssatutils))  # [dssatutils] was source(...)/' "$file"
  echo "rewired source() -> library(dssatutils) in $rel"
done

# 2) Guarded install + library at top of pipeline/config.R (only once).
CFG="pipeline/config.R"
if [ -f "$CFG" ] && ! grep -q "dssatutils" "$CFG"; then
  cp "$CFG" "$BACKUP/config.R"
  if [ "$USE_LOCAL" = "1" ]; then
    INSTALL_LINE="  remotes::install_local(\"$LOCAL_PKG\", force = FALSE, upgrade = \"never\")"
  else
    INSTALL_LINE="  remotes::install_github(\"$PKG_TAG\")"
  fi
  cat > /tmp/dssatutils_loader.R <<EOF
# --- Shared weather/soil utilities, extracted to the dssatutils package --------
# (Replaces the per-script source() of r_scripts/weather_*.R & soil_*.R.)
if (!requireNamespace("dssatutils", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
$INSTALL_LINE
}
suppressMessages(library(dssatutils))
# ------------------------------------------------------------------------------

EOF
  cat /tmp/dssatutils_loader.R "$CFG" > /tmp/config_new.R && mv /tmp/config_new.R "$CFG"
  echo "prepended dssatutils loader to $CFG"
else
  echo "config.R missing or already patched — skipping loader insert"
fi

# 3) Delete the now-duplicated weather/soil scripts (KEEP landcover_*).
for f in r_scripts/weather_daymet.R r_scripts/weather_gridmet.R \
         r_scripts/weather_nasapower.R r_scripts/soil_soilgrids.R \
         r_scripts/soil_soilgrids_online.R r_scripts/soil_ssurgo.R; do
  if [ -f "$f" ]; then
    mkdir -p "$BACKUP/$(dirname "$f")"; cp "$f" "$BACKUP/$f"
    git rm -q "$f" 2>/dev/null || rm -f "$f"
    echo "removed $f"
  fi
done
echo "kept: r_scripts/landcover_raster.R, r_scripts/landcover_raster_to_gridpoints.R"

echo
echo "PHASE 2 done. Review:  git -C \"$REPO\" diff"
echo "Then in R (repo has no lockfile):  renv::init(); renv::snapshot()"
echo "Backups: $BACKUP"

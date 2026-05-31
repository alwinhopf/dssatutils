#!/usr/bin/env bash
# Smoke-test that dssatutils installs and imports locally, BEFORE relying on the
# GitHub remote. Safe to run repeatedly. Writes a report to /tmp/dssatutils_check.txt.
set -uo pipefail
cd "$(dirname "$0")/.."
PKG="$(pwd)"
OUT=/tmp/dssatutils_check.txt
: > "$OUT"

echo "== Python: import via sys.path (no install) ==" | tee -a "$OUT"
python3 -c "import sys; sys.path.insert(0,'$PKG/python'); import dssatutils; print('import OK; n=', len(dssatutils.__all__))" 2>&1 | tee -a "$OUT"

echo "== Python: editable install (optional, comment out if undesired) ==" | tee -a "$OUT"
# python3 -m pip install -e "$PKG" 2>&1 | tail -5 | tee -a "$OUT"

echo "== Python: each submodule imports (deps permitting) ==" | tee -a "$OUT"
for m in weather_daymet weather_gridmet weather_nasapower weather_agera5 \
         weather_openmeteo weather_nasapower_chirps soil_soilgrids \
         soil_soilgrids_online soil_ssurgo soil_hwsd; do
  python3 -c "import sys; sys.path.insert(0,'$PKG/python'); import importlib; importlib.import_module('dssatutils.$m'); print('  ok: $m')" 2>&1 | tee -a "$OUT"
done

echo "== R: source each file + confirm exported fn exists ==" | tee -a "$OUT"
Rscript -e '
  files <- list.files("R", pattern="\\.R$", full.names=TRUE)
  for (f in files) { ok <- tryCatch({ sys.source(f, envir=new.env()); "ok" },
                                     error=function(e) paste("ERR:", conditionMessage(e)))
                     cat(sprintf("  %-32s %s\n", basename(f), ok)) }
' 2>&1 | tee -a "$OUT"

echo "Report written to $OUT"

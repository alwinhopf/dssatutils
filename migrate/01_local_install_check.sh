#!/usr/bin/env bash
# Smoke-test that dssatutils is syntactically sound and (deps permitting) imports.
# Safe to run repeatedly. Writes a report to /tmp/dssatutils_check.txt.
#
# Two tiers:
#   TIER 1 (always works, no heavy deps): byte-compile every .py and parse every .R
#   TIER 2 (optional, needs numpy/pandas/geo stack): real import of each submodule
set -uo pipefail
cd "$(dirname "$0")/.."
PKG="$(pwd)"
OUT=/tmp/dssatutils_check.txt
: > "$OUT"

say() { echo "$@" | tee -a "$OUT"; }

say "== TIER 1a: Python byte-compile (syntax only, no imports) =="
if python3 -m py_compile "$PKG"/python/dssatutils/*.py 2>>"$OUT"; then
  say "  OK: all Python modules compile"
else
  say "  FAIL: a Python module has a syntax error (see above)"
fi

say "== TIER 1b: dssatutils package object imports (lazy, cheap) =="
python3 -c "import sys; sys.path.insert(0,'$PKG/python'); import dssatutils; assert len(dssatutils.__all__)==11; print('  OK: import dssatutils, __all__ has', len(dssatutils.__all__), 'fns')" 2>>"$OUT" \
  || say "  FAIL: could not import dssatutils package object"

say "== TIER 1c: R parse every file (syntax only, no library() exec) =="
Rscript -e 'ok<-TRUE; for (f in list.files("R", pattern="[.]R$", full.names=TRUE)) { r<-tryCatch({parse(file=f); "ok"}, error=function(e) paste("ERR:",conditionMessage(e))); cat(sprintf("  %-30s %s\n", basename(f), r)); if(r!="ok") ok<-FALSE }; if(!ok) quit(status=1)' 2>>"$OUT" \
  || say "  FAIL: an R file failed to parse"

say "== TIER 2 (optional): real submodule imports — needs full geo stack =="
say "   Skipped by default. To run: SET DEEP=1 ./migrate/01_local_install_check.sh"
if [ "${DEEP:-0}" = "1" ]; then
  for m in weather_daymet weather_gridmet weather_nasapower weather_agera5 \
           weather_openmeteo weather_nasapower_chirps soil_soilgrids \
           soil_soilgrids_online soil_ssurgo soil_ssurgo_alderman soil_hwsd; do
    python3 -c "import sys; sys.path.insert(0,'$PKG/python'); import importlib; importlib.import_module('dssatutils.$m'); print('  ok: $m')" 2>>"$OUT" \
      || say "  (import failed: $m — likely a missing optional dep in this env)"
  done
fi

say "Report written to $OUT"

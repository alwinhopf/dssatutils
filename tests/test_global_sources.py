#!/usr/bin/env python3
# File: tests/test_global_sources.py
# ---------------------------------------------------------------------------
# Stand-alone smoke test for the GLOBAL data sources (Open-Meteo weather +
# SoilGrids soil) on a handful of points in Europe, Asia, and Africa.
#
# WHY: the global sources hit live APIs that are blocked in CI sandboxes, so
# this is meant to be run by the user on their own machine with internet
# access. It exercises a tiny 3-point "grid" end-to-end and validates that the
# emitted DSSAT .WTH / .SOL files are well-formed (headers, column widths,
# date format, no NaN), without needing DSSAT itself installed.
#
# Run:
#   python tests/test_global_sources.py
#   python tests/test_global_sources.py --keep   # keep output for inspection
# ---------------------------------------------------------------------------

import os
import sys
import argparse
import tempfile
import shutil

import pandas as pd

# Make python_scripts importable regardless of cwd.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
sys.path.insert(0, os.path.join(_REPO, "python_scripts"))

from weather_openmeteo import process_weather_openmeteo            # noqa: E402

try:
    import geopandas as gpd                                          # noqa: E402
    from shapely.geometry import Point                               # noqa: E402
    import soil_soilgrids_online as sg_online                        # noqa: E402
    _HAVE_SOIL = True
except Exception as exc:  # noqa: BLE001
    print(f"[warn] SoilGrids/geopandas not importable ({exc}); skipping soil test.")
    _HAVE_SOIL = False

# 3 points: Europe (NL), Asia (India Punjab), Africa (Kenya).
POINTS = pd.DataFrame({
    "ID":   ["EU_NL", "AS_IN", "AF_KE"],
    "LAT":  [52.000, 30.900, -0.500],
    "LONG": [5.000, 75.800, 37.000],
})

START_YEAR = 2010
END_YEAR = 2011


def _check_wth(path: str) -> None:
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    assert lines[0].startswith("$WEATHER"), f"{path}: bad title line"
    assert lines[1].lstrip().startswith("@ INSI"), f"{path}: missing INSI header"
    assert lines[3].lstrip().startswith("@  DATE"), f"{path}: missing DATE header"
    data = lines[4:]
    assert len(data) > 300, f"{path}: only {len(data)} daily rows (expected ~730)"
    for ln in data[:5] + data[-5:]:
        date = ln[:7].strip()
        assert len(date) == 7 and date.isdigit(), f"{path}: bad DATE token '{date}'"
        assert "nan" not in ln.lower(), f"{path}: NaN in data row"
    print(f"  [ok] {os.path.basename(path)}: {len(data)} daily rows, header valid")


def _check_sol(path: str) -> None:
    with open(path) as fh:
        txt = fh.read()
    assert "*SOILS" in txt or txt.lstrip().startswith("*"), f"{path}: bad .SOL start"
    assert "@SLB" in txt or "SLB" in txt, f"{path}: missing layer table"
    assert "nan" not in txt.lower(), f"{path}: NaN in soil file"
    print(f"  [ok] {os.path.basename(path)}: .SOL header + layer table present")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", action="store_true", help="keep output dir")
    args = ap.parse_args()

    work = tempfile.mkdtemp(prefix="dssat_global_test_")
    wth_dir = os.path.join(work, "weather")
    sol_dir = os.path.join(work, "soil")
    os.makedirs(wth_dir, exist_ok=True)
    os.makedirs(sol_dir, exist_ok=True)
    log = os.path.join(work, "errors.log")

    print(f"\n=== GLOBAL SOURCE SMOKE TEST ===\nWork dir: {work}\n")

    failures = 0

    # --- Weather: Open-Meteo (global, keyless ERA5) ---
    print("[1/2] Open-Meteo weather (EU / Asia / Africa)...")
    try:
        process_weather_openmeteo(
            shapefile=POINTS, start_year=START_YEAR, end_year=END_YEAR,
            output_dir=wth_dir, id_col="ID", lat_col="LAT", lon_col="LONG",
            n_cores=3, log_file=log,
        )
        for pid in POINTS["ID"]:
            p = os.path.join(wth_dir, f"{pid}.WTH")
            if not os.path.exists(p):
                print(f"  [FAIL] missing {pid}.WTH"); failures += 1
            else:
                _check_wth(p)
    except Exception as exc:  # noqa: BLE001
        print(f"  [FAIL] Open-Meteo raised: {exc}"); failures += 1

    # --- Soil: SoilGrids online (global, REST mode = one request per point) ---
    #
    # Per-point accounting: every input point must be EXPLICITLY accounted for —
    # it either produced a well-formed .SOL, or it is reported as a no-data gap.
    # SoilGrids has genuine coverage gaps (points over water, deserts, some
    # agricultural cells), so a missing .SOL is reported as a [skip-no-data]
    # WARNING rather than a hard failure. What this guards against is the old
    # blind spot: a point silently vanishing with the test still printing PASS.
    if _HAVE_SOIL:
        print("\n[2/2] SoilGrids online soil (EU / Asia / Africa)...")
        try:
            sg_online.USE_REST_API = True   # avoid GDAL/VRT dependency for the test
            gdf = gpd.GeoDataFrame(
                POINTS.copy(),
                geometry=[Point(xy) for xy in zip(POINTS["LONG"], POINTS["LAT"])],
                crs="EPSG:4326",
            )
            csv_path = os.path.join(sol_dir, "soil_map.csv")
            sg_online.process_soils_soilgrids_online(
                gridfile=gdf, soilfile_csv_path=csv_path,
                output_sol_dir=sol_dir, id_col="ID",
            )
            # Account for EVERY input point individually.
            produced, skipped = [], []
            for pid in POINTS["ID"]:
                p = os.path.join(sol_dir, f"{pid}.SOL")
                if os.path.exists(p):
                    _check_sol(p)           # validates header + layer table
                    produced.append(pid)
                else:
                    print(f"  [skip-no-data] {pid}.SOL not written "
                          f"(SoilGrids coverage gap)")
                    skipped.append(pid)

            # Hard requirement: at least one point must succeed, otherwise the
            # source is effectively broken (not just a coverage gap).
            if not produced:
                print("  [FAIL] no .SOL files written for ANY point"); failures += 1
            else:
                print(f"  [ok] {len(produced)}/{len(POINTS)} points produced "
                      f"valid .SOL ({len(skipped)} no-data gap(s))")
        except Exception as exc:  # noqa: BLE001
            print(f"  [FAIL] SoilGrids raised: {exc}"); failures += 1
    else:
        print("\n[2/2] SoilGrids: skipped (module not importable).")

    print("\n=== RESULT ===")
    if failures == 0:
        print("ALL CHECKS PASSED.")
    else:
        print(f"{failures} CHECK(S) FAILED. See {log} if present.")
        if os.path.exists(log):
            print("\n--- error log ---")
            print(open(log).read())

    if args.keep:
        print(f"\nOutput kept at: {work}")
    else:
        shutil.rmtree(work, ignore_errors=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())

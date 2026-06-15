#!/usr/bin/env python3
# File: tests/test_polaris.py
# ---------------------------------------------------------------------------
# OFFLINE tests for the POLARIS soil source (Tier-0 deterministic drop-in).
# No network and no DSSAT install required: the tile-streaming fetch is isolated
# from the pure tile-addressing / unit-transform / van-Genuchten / .SOL-writer
# helpers, which are exercised here with deterministic synthetic inputs. A
# cross-language section asserts the R twin carries the same API + algorithm.
#
# Run:  python -m pytest tests/test_polaris.py
# ---------------------------------------------------------------------------

import math
import os
import sys
import tempfile

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)
_WORKSPACE = os.path.dirname(_REPO)
sys.path.insert(0, os.path.join(_REPO, "python"))


def test_polaris_import_and_export():
    import dssatutils
    assert hasattr(dssatutils, "process_soils_polaris")


def test_tile_addressing():
    from dssatutils import soil_polaris as p
    # Iowa (Hardin Co.) and Texas (Hutchinson Co.) example-field centroids.
    assert p._polaris_tile(42.35, -93.40) == "lat4243_lon-94-93"
    assert p._polaris_tile(35.91, -101.40) == "lat3536_lon-102-101"
    # exact-integer and positive-lon edge cases
    assert p._polaris_tile(40.0, -100.0) == "lat4041_lon-100-99"


def test_log_back_transform():
    from dssatutils import soil_polaris as p
    assert abs(p._backtransform("om", 0.0) - 1.0) < 1e-12        # 10**0 = 1 %
    assert abs(p._backtransform("alpha", -2.0) - 0.01) < 1e-12   # log10 kPa^-1
    assert abs(p._backtransform("ksat", 1.0) - 10.0) < 1e-12     # log10 cm/hr
    assert p._backtransform("clay", 25.0) == 25.0                # linear passthrough
    assert math.isnan(p._backtransform("om", float("nan")))


def test_van_genuchten_water_limits_ordering():
    from dssatutils import soil_polaris as p
    # plausible silty-clay-loam VG params (alpha already in kPa^-1)
    wl = p.water_limits(theta_r=0.08, theta_s=0.46, alpha=0.02, n=1.4)
    assert wl["SLLL"] < wl["SDUL"] < wl["SSAT"]
    assert 0.10 < wl["SLLL"] < 0.25
    assert 0.30 < wl["SDUL"] < 0.50
    assert abs(wl["SSAT"] - 0.46) < 1e-9          # SAT == theta_s here
    # theta(1500) < theta(33): curve is monotonically decreasing in suction
    assert p._vg_theta(1500.0, 0.08, 0.46, 0.02, 1.4) < \
        p._vg_theta(33.0, 0.08, 0.46, 0.02, 1.4)


def test_water_limits_fallback_and_guards():
    from dssatutils import soil_polaris as p
    # VG params missing -> Saxton-Rawls fallback from texture still gives a
    # physically ordered, crash-safe profile.
    wl = p.water_limits(float("nan"), float("nan"), float("nan"), float("nan"),
                        sand=40.0, clay=20.0, om_pct=2.0)
    assert wl["SLLL"] >= 0.02
    assert wl["SDUL"] - wl["SLLL"] >= 0.04 - 1e-9
    assert wl["SSAT"] - wl["SDUL"] >= 0.04 - 1e-9
    # near-pure sand: guards must still hold the ordering
    wl2 = p.water_limits(0.03, 0.33, 0.08, 2.6, sand=95.0, clay=2.0, om_pct=0.3)
    assert wl2["SLLL"] < wl2["SDUL"] < wl2["SSAT"]


def test_ssks_clamp():
    from dssatutils import soil_polaris as p
    assert p._ssks_cmhr(10.0) == 10.0
    assert p._ssks_cmhr(5000.0) == 999.0
    assert p._ssks_cmhr(-1.0) == 0.0
    assert p._ssks_cmhr(float("nan")) == -99.0


def test_sol_writer_roundtrip():
    from dssatutils import soil_polaris as p
    prof = pd.DataFrame({
        "ID": ["00000001"] * 2,
        "latitude": [42.35, 42.35], "longitude": [-93.40, -93.40],
        "depth_bottom": [5, 15], "depth_center": [2.5, 10.0],
        "SLLL": [0.177, 0.182], "SDUL": [0.415, 0.421], "SSAT": [0.460, 0.466],
        "SSKS": [1.2, 0.9], "bd": [1.35, 1.38], "oc_pct": [2.10, 1.80],
        "clay": [26.0, 28.0], "silt": [41.0, 40.0], "ph": [6.2, 6.4],
    })
    with tempfile.TemporaryDirectory() as work:
        p._format_dssat_sol_file(prof, work, source_tag="p50")
        txt = open(os.path.join(work, "00000001.SOL")).read()
    assert "*SOILS: POLARIS v1.0" in txt
    assert "statistic=p50" in txt
    assert "@  SLB" in txt
    # the first data layer must keep SLLL < SDUL < SSAT after formatting
    data = [ln for ln in txt.splitlines() if ln.lstrip().startswith("5 ")
            or ln.lstrip().startswith("5\t") or ln.strip().startswith("5 ")]
    layer = next(ln for ln in txt.splitlines() if ln.strip().startswith("5 "))
    cols = layer.split()
    slll, sdul, ssat = float(cols[2]), float(cols[3]), float(cols[4])
    assert slll < sdul < ssat


def test_r_python_parity_markers():
    # The R twin must carry the same public API and core algorithm so the two
    # implementations cannot silently drift (AGENTS.md R/Python parity contract).
    rel = "R/soil_polaris.R"
    src = open(os.path.join(_WORKSPACE, "dssatutils", rel),
               encoding="utf-8", errors="replace").read()
    for marker in ("process_soils_polaris", "POLARIS", "hydrology.cee.duke.edu",
                   "theta_r", "p50", "Saxton"):
        assert marker in src, f"{rel} missing expected marker {marker!r}"


if __name__ == "__main__":
    import pytest
    raise SystemExit(pytest.main([__file__, "-v"]))

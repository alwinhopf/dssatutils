"""Offline regressions read the same columns as DSSAT PARSE_HEADERS."""
from pathlib import Path
import importlib
import re

import numpy as np
import pandas as pd
import pytest

from dssatutils.weather_format import format_wth_value, wind_run
from dssatutils.weather_repair import _parse_daily_rows, _parse_wth, _write_wth


def _native_values(header, line):
    tokens = list(re.finditer(r"\S+", header))
    # DSSAT starts each field one column beyond the preceding header's end.
    return {m.group(): float(line[tokens[i-1].end()+1:m.end()])
            for i, m in enumerate(tokens) if i and m.group() != "DATE"}


def test_numeric_fields_preserve_native_separator_and_rounding():
    fixture = pd.read_csv(Path(__file__).parent / "fixtures/dssat_numeric_fields.csv", keep_default_na=False, dtype={"expected": str})
    for row in fixture.itertuples():
        field = format_wth_value(row.value, row.decimals)
        assert len(field) == 6 and field.startswith(" ")
        assert field.strip() == row.expected
    with pytest.raises(ValueError, match="cannot fit"):
        format_wth_value(100000)
    assert [wind_run(v) for v in (3, -99, np.nan, np.inf)] == pytest.approx([259.2, -99, -99, -99])


@pytest.mark.parametrize("provider", ["cmfd", "dwd", "eobs", "xavier", "gridded_common", "agera5", "era5land"])
def test_provider_wind_units_and_native_columns(provider, tmp_path):
    module = importlib.import_module("dssatutils.weather_" + provider)
    frame = pd.DataFrame(dict(DATE=["2024001", "2024002"], YEAR=2024, MM=1, DOY=[1, 2],
                              SRAD=12., TMAX=-5., TMIN=-12.3, RAIN=[999.96, 1200.],
                              TDEW=-13., RH2M=76.9, WIND=[13., -99.]))
    if provider == "dwd":
        module._write_wth(frame, "TEST", 50., 10., 100., str(tmp_path))
    elif provider == "gridded_common":
        module.write_wth(frame, "TEST", 50., 10., str(tmp_path), "test", "TEST")
    elif provider == "era5land":
        module._write_dssat_weather_file(frame, 50., 10., str(tmp_path / "TEST.WTH"), "TEST")
    else:
        module._write_wth(frame, "TEST", 50., 10., str(tmp_path))
    lines = (tmp_path / "TEST.WTH").read_text().splitlines()
    header = next(l for l in lines if l.startswith("@  DATE"))
    rows = [l for l in lines if l[:7].isdigit()]
    values = [_native_values(header, l) for l in rows]
    assert len(rows) == 2 and all(len(l) == 49 for l in rows)
    assert [v["WIND"] for v in values] == [1123., -99.]
    assert [v["RAIN"] for v in values] == [1000., 1200.]
    assert all(v["TMIN"] == -12.3 and v["RH2M"] == 76.9 for v in values)
    if provider == "era5land":
        assert lines[2].split()[-1] == "10.0"


@pytest.mark.parametrize("provider", ["isdasoil", "lucas", "ssurgo", "gnatsgo"])
def test_soil_conductivity_rounding_and_native_columns(provider, tmp_path):
    module = importlib.import_module("dssatutils.soil_" + provider)
    frame = pd.DataFrame(dict(ID="TEST", latitude=50., longitude=10., depth_bottom=[20, 50],
                              SLLL=.1, SDUL=.2, SSAT=.4, SSKS=["99.996", "150"],
                              bulk_density=1.4, om_pct=2., clay_pct=20., silt_pct=30.))
    module._write_sol(frame, str(tmp_path))
    lines = (tmp_path / "TEST.SOL").read_text().splitlines()
    i = next(i for i,l in enumerate(lines) if l.startswith("@  SLB"))
    values = [_native_values(lines[i], l) for l in lines[i+1:] if l.strip()]
    assert [v["SSKS"] for v in values] == [100., 150.]
    assert [v["SRGF"] for v in values] == [1., 1.]


def test_legacy_adjacent_rows_survive_repair_for_native_reader(tmp_path):
    path = tmp_path / "TEST.WTH"
    path.write_bytes((Path(__file__).parent / "fixtures/weather_adjacent_wind.txt").read_bytes())
    lines, index, _, data, _ = _parse_wth(path)
    assert len(data) == 7 and data.WIND.eq(1062.7).all()
    _write_wth(path, lines, index, data)
    output = path.read_text().splitlines()
    assert all(_native_values(output[index], l)["WIND"] == 1063 for l in output[index+1:])
    with pytest.raises(ValueError, match="Malformed"):
        _parse_daily_rows(["2024001 12 broken"])


def test_alderman_conductivity_does_not_shift_following_fields(tmp_path):
    from dssatutils.soil_ssurgo_alderman import _write_dssat_soil_file
    layers = pd.DataFrame(dict(SLB=[20, 50], SLMH="-99", SLLL=.1, SDUL=.2, SSAT=.4,
                               SRGF=1., SSKS=["99.996", "150"], SBDM=1.4, SLOC=1.2,
                               SLCL=20., SLSI=30., SLCF=0., SLNI=-99., SLHW=6.,
                               SLHB=-99., SCEC=-99., SADC=-99.))
    profile = dict(profile_id="TEST", site="TEST", country="USA", latitude=50., longitude=10.,
                   scs_family="-99", scom="BN", salb=.13, slu1=6., sldr=.6, slro=73.,
                   slnf=1., slpf=1., smhb="IB001", smpx="IB001", smke="IB001", layers=layers)
    _write_dssat_soil_file(profile, str(tmp_path))
    lines = (tmp_path / "TEST.SOL").read_text().splitlines()
    i = next(i for i,l in enumerate(lines) if l.startswith("@  SLB"))
    values = [_native_values(lines[i], l) for l in lines[i+1:] if l.strip()]
    assert [v["SSKS"] for v in values] == [100., 150.]
    assert [v["SBDM"] for v in values] == [1.4, 1.4]

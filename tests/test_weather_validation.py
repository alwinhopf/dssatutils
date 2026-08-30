from pathlib import Path

import pandas as pd

from dssatutils import is_wth_valid
from dssatutils.weather_agera5 import _write_wth
from dssatutils.weather_nasapower import _NASA_PARAMS, _fetch_nasa_power


def _write_sample(path: Path, rows: list[str]) -> None:
    path.write_text(
        "$WEATHER DATA: test\n"
        "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT\n"
        "  TEST   0.0000   0.0000   -99  10.0  20.0   2.0  10.0\n"
        "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND\n" +
        "\n".join(rows) + "\n",
        encoding="utf-8",
    )


def test_fixed_width_adjacent_negative_values_are_valid(tmp_path):
    path = tmp_path / "TEST.WTH"
    rows = [
        f"{'2024001':>7}{12.0:6.1f}{-10.2:6.1f}{-12.3:6.1f}{0.0:6.1f}{-15.0:6.1f}{40.0:6.1f}{3.0:6.1f}",
        f"{'2024002':>7}{13.0:6.1f}{-9.0:6.1f}{-11.0:6.1f}{0.0:6.1f}{-14.0:6.1f}{42.0:6.1f}{3.2:6.1f}",
    ]
    _write_sample(path, rows)
    assert is_wth_valid(path, end_year=2024)


def test_weather_validator_rejects_date_gaps(tmp_path):
    path = tmp_path / "TEST.WTH"
    _write_sample(path, [
        "2024001 12.0 10.0 1.0 0.0 0.0 40.0 3.0",
        "2024003 12.0 10.0 1.0 0.0 0.0 40.0 3.0",
    ])
    assert not is_wth_valid(path, end_year=2024)


def test_weather_validator_rejects_absolute_zero_temperature(tmp_path):
    path = tmp_path / "TEST.WTH"
    _write_sample(path, [
        f"{'2018001':>7}{0.0:6.1f}{-273.1:6.1f}{-273.1:6.1f}{0.0:6.1f}{-273.1:6.1f}{0.0:6.1f}{0.0:6.1f}"
    ])
    assert not is_wth_valid(path, end_year=2018)


def test_weather_validator_can_require_complete_core_forcing(tmp_path):
    path = tmp_path / "TEST.WTH"
    _write_sample(path, [
        "2024001 -99 10.0 1.0 0.0 0.0 40.0 -99",
        "2024002 12.0 10.0 1.0 0.0 0.0 40.0 -99",
    ])
    assert is_wth_valid(path, end_year=2024)
    assert not is_wth_valid(
        path, end_year=2024, required_columns=("SRAD", "TMAX", "TMIN", "RAIN")
    )


def test_weather_validator_can_require_all_agera5_forcing(tmp_path):
    path = tmp_path / "TEST.WTH"
    _write_sample(path, [
        "2024001 12.0 10.0 1.0 0.0 -99 -99 -99",
        "2024002 12.0 10.0 1.0 0.0 -99 -99 -99",
    ])
    core = ("SRAD", "TMAX", "TMIN", "RAIN")
    agera5 = core + ("TDEW", "RH2M", "WIND")

    assert is_wth_valid(path, end_year=2024, required_columns=core)
    assert not is_wth_valid(path, end_year=2024, required_columns=agera5)


def test_agera5_writer_defers_physical_validation_to_shared_validator(tmp_path):
    frame = pd.DataFrame({
        "DATE": ["2018001", "2018002"], "YEAR": [2018, 2018], "MM": [1, 1],
        "SRAD": [12.0, 12.0], "TMAX": [5.8, 10.0],
        "TMIN": [6.0, 2.0], "RAIN": [0.0, 0.0],
        "TDEW": [0.0, 0.0], "RH2M": [60.0, 60.0], "WIND": [3.0, 3.0],
    })
    path = Path(_write_wth(frame, "TEST", 33.7, -102.5, str(tmp_path)))

    assert path.exists()
    assert "   5.8   6.0" in path.read_text(encoding="utf-8")
    assert not is_wth_valid(
        path,
        end_year=2018,
        required_columns=("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"),
    )


def test_nasa_power_normalizes_json_null_to_numeric_missing(monkeypatch):
    class Response:
        def raise_for_status(self):
            return None

        def json(self):
            parameters = {
                name: {"19901103": (None if name == "ALLSKY_SFC_SW_DWN" else 1.0)}
                for name in _NASA_PARAMS
            }
            return {"properties": {"parameter": parameters}}

    monkeypatch.setattr(
        "dssatutils.weather_nasapower.requests.get",
        lambda *args, **kwargs: Response(),
    )

    frame = _fetch_nasa_power(34.57, -102.60, "19901103", "19901103", retries=1)

    assert frame.loc[0, "ALLSKY_SFC_SW_DWN"] == -99.0
    assert frame[_NASA_PARAMS].notna().all().all()

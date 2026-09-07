import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))

from dssatutils.weather_repair import (  # noqa: E402
    audit_weather_file_quality,
    repair_weather_file_date_gaps,
    repair_weather_file_temperature_inversions,
)


def _write_wth(path: Path) -> None:
    path.write_text(
        "\n".join([
            "$WEATHER DATA: TEST",
            "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
            "  TEST   0.0000   0.0000   -99  20.0  10.0   2.0   2.0",
            "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
            "2024001  15.0  20.0  10.0   0.0   8.0  60.0   2.0",
            "2024002  15.0  22.0  12.0   0.0   9.0  60.0   2.0",
            "2024003  15.0   5.0  15.0   0.0  10.0  60.0   2.0",
            "2024004  15.0  24.0  14.0   0.0  11.0  60.0   2.0",
            "2024005  15.0  26.0  16.0   0.0  12.0  60.0   2.0",
            "2024006  15.0  28.0  18.0   0.0  13.0  60.0   2.0",
        ]) + "\n",
        encoding="utf-8",
    )


def test_repair_weather_temperature_inversion_uses_neighbor_means(tmp_path):
    wth = tmp_path / "00000001.WTH"
    log = tmp_path / "weather_repair.log"
    _write_wth(wth)

    summary = repair_weather_file_temperature_inversions(
        wth,
        max_gap_days=3,
        window_days=2,
        log_file=log,
    )

    assert summary.loc[0, "status"] == "repaired"
    assert int(summary.loc[0, "repaired_count"]) == 1

    rows = []
    for line in wth.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("2024"):
            rows.append(line.split())
    dat = pd.DataFrame(rows, columns=["DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"])
    dat[["TMAX", "TMIN"]] = dat[["TMAX", "TMIN"]].astype(float)

    repaired = dat.loc[dat["DATE"] == "2024003"].iloc[0]
    assert repaired["TMAX"] == 23.0
    assert repaired["TMIN"] == 13.0
    assert repaired["TMIN"] <= repaired["TMAX"]
    assert "issue=TMIN_GT_TMAX status=repaired" in log.read_text(encoding="utf-8")


def test_repair_weather_temperature_inversion_swap(tmp_path):
    wth = tmp_path / "00000001.WTH"
    log = tmp_path / "weather_repair.log"
    # Create file with one small inversion (diff <= 2.0) and one large inversion (diff > 2.0)
    wth.write_text(
        "\n".join([
            "$WEATHER DATA: TEST",
            "@ INSI      LAT     LONG  ELEV   TAV   AMP REFHT WNDHT",
            "  TEST   0.0000   0.0000   -99  20.0  10.0   2.0   2.0",
            "@  DATE  SRAD  TMAX  TMIN  RAIN  TDEW  RH2M  WIND",
            "2024001  15.0  20.0  10.0   0.0   8.0  60.0   2.0",
            "2024002  15.0  14.0  15.0   0.0   9.0  60.0   2.0",  # diff = 1.0 <= 2.0 -> swap
            "2024003  15.0  10.0  15.0   0.0  10.0  60.0   2.0",  # diff = 5.0 > 2.0 -> unrepaired
            "2024004  15.0  24.0  14.0   0.0  11.0  60.0   2.0",
        ]) + "\n",
        encoding="utf-8",
    )

    summary = repair_weather_file_temperature_inversions(
        wth,
        method="swap",
        max_inversion_c=2.0,
        log_file=log,
    )

    assert summary.loc[0, "status"] == "repaired"
    assert int(summary.loc[0, "repaired_count"]) == 1
    assert int(summary.loc[0, "unrepaired_count"]) == 1

    rows = []
    for line in wth.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("2024"):
            rows.append(line.split())
    dat = pd.DataFrame(rows, columns=["DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"])
    dat[["TMAX", "TMIN"]] = dat[["TMAX", "TMIN"]].astype(float)

    day2 = dat.loc[dat["DATE"] == "2024002"].iloc[0]
    assert day2["TMAX"] == 15.0
    assert day2["TMIN"] == 14.0

    day3 = dat.loc[dat["DATE"] == "2024003"].iloc[0]
    assert day3["TMAX"] == 10.0
    assert day3["TMIN"] == 15.0



def test_repair_weather_date_gap_inserts_neighbor_mean_row(tmp_path):
    wth = tmp_path / "00000002.WTH"
    log = tmp_path / "weather_repair.log"
    _write_wth(wth)
    lines = [line for line in wth.read_text(encoding="utf-8").splitlines() if not line.startswith("2024003")]
    wth.write_text("\n".join(lines) + "\n", encoding="utf-8")

    summary = repair_weather_file_date_gaps(
        wth,
        max_gap_days=3,
        window_days=2,
        log_file=log,
    )

    assert summary.loc[0, "status"] == "repaired"
    assert int(summary.loc[0, "repaired_count"]) == 1
    rows = []
    for line in wth.read_text(encoding="utf-8").splitlines():
        if line.strip().startswith("2024"):
            rows.append(line.split())
    dat = pd.DataFrame(rows, columns=["DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"])
    dat[["TMAX", "TMIN"]] = dat[["TMAX", "TMIN"]].astype(float)
    inserted = dat.loc[dat["DATE"] == "2024003"].iloc[0]
    assert inserted["TMAX"] == 23.0
    assert inserted["TMIN"] == 13.0
    assert "issue=DATE_GAP status=repaired" in log.read_text(encoding="utf-8")


def test_audit_weather_quality_flags_suspicious_rows(tmp_path):
    wth = tmp_path / "00000003.WTH"
    _write_wth(wth)

    audit = audit_weather_file_quality(wth, flatline_days=3)

    assert "tmin_gt_tmax" in set(audit["issue"])
    assert "RAIN_flatline" in set(audit["issue"])

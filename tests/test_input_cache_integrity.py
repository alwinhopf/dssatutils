"""Offline regressions for cached SOL formatting and GRIDMET point identity."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from dssatutils import rebuild_soil_files_from_mapping, soil_file_issue
from dssatutils.weather_gridmet import _gridmet_extract_cells

FIXTURES = Path(__file__).parent / "fixtures"


def test_gridmet_retains_invalid_duplicate_and_out_of_order_points():
    points = pd.read_csv(FIXTURES / "gridmet_point_alignment.csv")
    cells = points.cell.fillna(0).astype(int).to_numpy() - 1
    data = np.array([[[10, 20, 30]], [[11, 21, 31]]])
    result = _gridmet_extract_cells(data, np.where(cells < 0, -1, 0), cells, 2)
    np.testing.assert_allclose(result, points[["day1", "day2"]].to_numpy(), equal_nan=True)
    assert np.isnan(_gridmet_extract_cells(data, [-1], [-1], 2)).all()
    np.testing.assert_array_equal(_gridmet_extract_cells(data, [0], [1], 1), [[20]])


@pytest.mark.parametrize("source", ["SSURGO", "GNATSGO"])
def test_soil_rebuild_fixed_columns_and_reject_shifted_cache(tmp_path, source):
    destination = tmp_path / source
    records = rebuild_soil_files_from_mapping(FIXTURES / "soil_mapping_rebuild.csv", destination, source)
    assert records.ID.tolist() == ["00000001"]
    path = Path(records.path.iloc[0])
    assert soil_file_issue(path) is None
    lines = path.read_text().splitlines()
    profile = lines[3]
    assert profile[1:11].strip() == "00000001"
    assert profile[13:24].strip().upper() == source
    assert float(profile[31:36]) == 200
    assert float(lines[5][25:33]) == 33.682
    assert float(lines[5][34:42]) == -89.478
    h = next(i for i, line in enumerate(lines) if line.startswith("@  SLB"))
    rows = lines[h + 1:h + 5]
    assert [int(row[:6]) for row in rows] == [5, 20, 100, 200]
    assert [float(row[13:18]) for row in rows] == [.1, .11, .12, .13]
    assert float(rows[0][37:42]) == 120.5
    assert float(rows[0][49:54]) == 1.5  # organic matter -> organic carbon, once
    original = path.read_bytes()
    with pytest.raises(ValueError, match="new directory"):
        rebuild_soil_files_from_mapping(FIXTURES / "soil_mapping_rebuild.csv", destination, source)
    assert path.read_bytes() == original
    for i in range(h + 2, h + 5):
        lines[i] = " " + lines[i]  # reproduce historic cat(vector) separators
    path.write_text("\n".join(lines))
    assert "fixed-width" in soil_file_issue(path)


def test_soil_rejects_missing_and_duplicate_depths(tmp_path):
    assert soil_file_issue(tmp_path / "missing.SOL") == "SOIL.SOL is missing"
    records = rebuild_soil_files_from_mapping(FIXTURES / "soil_mapping_rebuild.csv", tmp_path / "new", "SSURGO")
    path = Path(records.path.iloc[0])
    lines = path.read_text().splitlines()
    h = next(i for i, row in enumerate(lines) if row.startswith("@  SLB"))
    lines[h + 2] = lines[h + 1]
    path.write_text("\n".join(lines))
    assert "not strictly increasing" in soil_file_issue(path)


def test_alderman_full_width_horizon_conductivity_and_missing_values(tmp_path):
    from dssatutils.soil_ssurgo_alderman import _write_dssat_soil_file
    layers = pd.read_csv(FIXTURES / "soil_mapping_rebuild.csv").rename(columns={
        "bulk_density": "SBDM", "clay_pct": "SLCL", "silt_pct": "SLSI"})
    layers["SLB"] = [5, 20, 100, 200]
    layers["SLMH"] = "Ap/Bt"
    layers["SRGF"] = 1.
    layers["SSKS"] = [10.08, 120.5, 0.25, -99]
    layers["SLOC"] = layers.om_pct / 1.724
    for name in ["SLCF", "SLNI", "SLHW", "SLHB", "SCEC", "SADC"]:
        layers[name] = np.nan
    profile = dict(profile_id="00000001", site="00000001", country="USA", latitude=33.682,
                   longitude=-89.478, scs_family="", scom="SC", salb=.13, slu1=6., sldr=.6,
                   slro=73., slnf=1., slpf=1., smhb="IB001", smpx="IB001", smke="IB001", layers=layers)
    _write_dssat_soil_file(profile, str(tmp_path))
    path = tmp_path / "00000001.SOL"
    assert soil_file_issue(path) is None
    lines = path.read_text().splitlines()
    h = next(i for i, row in enumerate(lines) if row.startswith("@  SLB"))
    rows = lines[h + 1:h + 5]
    assert [row[7:12] for row in rows] == ["Ap/Bt"] * 4
    assert [float(row[37:42]) for row in rows] == [10.08, 120.5, .25, -99]
    for i in range(h + 1, h + 5):
        lines[i] = lines[i][1:]  # old Alderman writer used five-column SLB
    path.write_text("\n".join(lines))
    assert "fixed-width" in soil_file_issue(path)

"""Offline PRISM cache, completeness and native R/Python contracts."""
import copy
import json
from pathlib import Path
import subprocess
from unittest.mock import Mock

import numpy as np
import pandas as pd
import pytest

from dssatutils import weather_prism as prism
from dssatutils.discovery import find_rscript
from dssatutils.weather_validation import is_wth_valid

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def offline(monkeypatch):
    monkeypatch.setattr(prism.requests, "post", Mock(side_effect=AssertionError("Unexpected network access")))
    monkeypatch.setattr(prism.time, "sleep", lambda _: None)


def records(sdate="2020-01-01", edate="2020-12-31"):
    return [[d.strftime("%Y-%m-%d"), "26", "16", "T"] for d in pd.date_range(sdate, edate)]


def response(raw):
    result = Mock(status_code=200, text="fixture")
    result.json.return_value = {"data": raw}
    return result


def run_point(tmp_path, lat=30.5, lon=-84.5, pid="P1", **kwargs):
    pts = pd.DataFrame({"ID": [pid], "LAT": [lat], "LONG": [lon]})
    prism.process_weather_prism(pts, 2020, 2020, str(tmp_path / "wth"),
        "ID", "LAT", "LONG", 1, str(tmp_path / "prism.log"), str(tmp_path / "cache"), **kwargs)


def test_acis_parse_records_coercion_and_sentinels():
    raw = json.loads((ROOT / "tests/fixtures/prism_acis_records.json").read_text())["data"]
    df = prism._parse_acis_records(raw, 30.0)
    assert list(df["DATE"]) == ["2020001", "2020002", "2020005", "2020006"]
    assert list(df["RAIN"]) == [0.0, 5.4, -99.0, -99.0]
    assert df.iloc[1]["TMAX"] == 27.5
    assert np.isfinite(df["SRAD"]).all()
    assert (df["SRAD"] > 0).all()
    assert (prism._parse_acis_records(raw, 30.0, "none")["SRAD"] == -99).all()


def test_complete_cache_reused_offline_and_moved_point_downloads(tmp_path, monkeypatch):
    post = Mock(return_value=response(records()))
    monkeypatch.setattr(prism.requests, "post", post)
    run_point(tmp_path)
    assert post.call_count == 1
    assert is_wth_valid(tmp_path / "wth/P1.WTH", start_year=2020, end_year=2020,
                        required_columns=["TMAX", "TMIN", "RAIN", "SRAD"])
    post.side_effect = AssertionError("Cache should avoid network")
    run_point(tmp_path, pid="RENAMED")
    assert post.call_count == 1
    # Same ID at new coordinates must download, rather than relabel old forcing.
    post.side_effect = None
    moved = records()
    for row in moved:
        row[1] = "35"
    post.return_value = response(moved)
    run_point(tmp_path, lat=40.5, lon=-100.5)
    assert post.call_count == 2
    assert post.call_args.kwargs["json"]["loc"] == "-100.5000,40.5000"
    assert "35.0" in (tmp_path / "wth/P1.WTH").read_text()
    assert len(list((tmp_path / "cache/acis").glob("*.json"))) == 2
    assert not list((tmp_path / "cache/acis").glob(".prism-*"))


@pytest.mark.parametrize("mutation", ["short", "gap", "duplicate", "reordered", "missing_temp", "missing_rain", "nonfinite", "malformed", "empty"])
def test_incomplete_download_never_cached_or_published(tmp_path, monkeypatch, mutation):
    raw = records()
    if mutation == "short": raw = raw[:2]
    elif mutation == "gap": del raw[20]
    elif mutation == "duplicate": raw[20] = copy.deepcopy(raw[19])
    elif mutation == "reordered": raw[20], raw[21] = raw[21], raw[20]
    elif mutation == "missing_temp": raw[20][1] = "M"
    elif mutation == "missing_rain": raw[20][3] = "-999"
    elif mutation == "nonfinite": raw[20][1] = "Inf"
    elif mutation == "malformed": raw[20] = ["2020-01-21"]
    elif mutation == "empty": raw = []
    post = Mock(return_value=response(raw))
    monkeypatch.setattr(prism.requests, "post", post)
    run_point(tmp_path)
    assert post.call_count == 3
    assert not (tmp_path / "wth/P1.WTH").exists()
    assert not list((tmp_path / "cache").rglob("*.json"))
    assert "incomplete or invalid ACIS history" in (tmp_path / "prism.log").read_text()


@pytest.mark.parametrize("invalid", ["legacy", "coordinates", "units", "gap", "empty", "corrupt", "scalar"])
def test_invalid_cache_retried_and_replaced(tmp_path, monkeypatch, invalid):
    cache = tmp_path / "point.json"
    data = {"request": prism._acis_request(30.5, -84.5, "2020-01-01", "2020-01-03"),
            "data": records("2020-01-01", "2020-01-03")}
    if invalid == "legacy": del data["request"]
    elif invalid == "coordinates": data["request"]["loc"] = "-100.5000,40.5000"
    elif invalid == "units": data["request"]["units"] = "degreeF,degreeF,inch"
    elif invalid == "gap": del data["data"][1]
    elif invalid == "empty": data["data"] = []
    cache.write_text("{broken" if invalid == "corrupt" else json.dumps("invalid") if invalid == "scalar" else json.dumps(data))
    good = records("2020-01-01", "2020-01-03")
    post = Mock(return_value=response(good))
    monkeypatch.setattr(prism.requests, "post", post)
    assert prism._fetch_acis_point(30.5, -84.5, "2020-01-01", "2020-01-03", str(cache)) == good
    assert post.call_count == 1
    assert json.loads(cache.read_text())["request"] == prism._acis_request(30.5, -84.5, "2020-01-01", "2020-01-03")


def test_partial_response_retried_before_cache_write(tmp_path, monkeypatch):
    good = records("2020-01-01", "2020-01-03")
    post = Mock(side_effect=[response(good[:1]), response(good)])
    monkeypatch.setattr(prism.requests, "post", post)
    cache = tmp_path / "point.json"
    assert prism._fetch_acis_point(30.5, -84.5, "2020-01-01", "2020-01-03", str(cache)) == good
    assert post.call_count == 2
    assert json.loads(cache.read_text())["data"] == good


def test_existing_weather_preserved_on_failure(tmp_path, monkeypatch):
    path = tmp_path / "wth/P1.WTH"
    path.parent.mkdir()
    path.write_bytes(b"previous successful weather")
    monkeypatch.setattr(prism.requests, "post", Mock(return_value=response([])))
    run_point(tmp_path)
    assert path.read_bytes() == b"previous successful weather"


def test_assembler_independently_rejects_partial_history(tmp_path, monkeypatch):
    monkeypatch.setattr(prism, "_fetch_acis_point", lambda *a: records()[:2])
    run_point(tmp_path)
    assert not (tmp_path / "wth/P1.WTH").exists()


def test_parallel_points_publish_only_complete_histories(tmp_path, monkeypatch):
    monkeypatch.setattr(prism, "_fetch_acis_point", lambda lat, *a: records() if lat == 30.5 else records()[:2])
    points = pd.DataFrame({"ID": ["GOOD", "BAD"], "LAT": [30.5, 40.5], "LONG": [-84.5, -100.5]})
    prism.process_weather_prism(points, 2020, 2020, str(tmp_path / "wth"), "ID", "LAT", "LONG",
                                2, str(tmp_path / "prism.log"), str(tmp_path / "cache"))
    assert is_wth_valid(tmp_path / "wth/GOOD.WTH", start_year=2020, end_year=2020)
    assert not (tmp_path / "wth/BAD.WTH").exists()


def test_disabled_radiation_permits_only_optional_missing_forcing(tmp_path, monkeypatch):
    monkeypatch.setattr(prism.requests, "post", Mock(return_value=response(records())))
    run_point(tmp_path, srad_method="none")
    path = tmp_path / "wth/P1.WTH"
    assert is_wth_valid(path, start_year=2020, end_year=2020, required_columns=["TMAX", "TMIN", "RAIN"])
    assert not is_wth_valid(path, required_columns=["SRAD"])


def test_r_python_share_cache_and_wth_behavior(tmp_path, monkeypatch):
    rscript = find_rscript()
    if not rscript:
        pytest.skip("Rscript unavailable")
    probe = subprocess.run([rscript, "--vanilla", "-e",
        "quit(status=if(requireNamespace('jsonlite',quietly=TRUE)&&requireNamespace('httr',quietly=TRUE)) 0 else 1)"],
        capture_output=True, timeout=30)
    if probe.returncode:
        pytest.skip("Optional R JSON/HTTP dependencies unavailable")
    monkeypatch.setattr(prism.requests, "post", Mock(return_value=response(records())))
    run_point(tmp_path)
    # Read Python's cache in R, regenerate WTH, then read R's serialization in Python.
    script = tmp_path / "parity.R"
    script.write_text('''a <- commandArgs(TRUE)
for (name in c("weather_prism.R", "weather_gridded_common.R", "weather_format.R", "weather_solar.R", "weather_validation.R", "utils.R")) source(file.path(a[1], name))
assignInNamespace("POST", function(...) stop("Network must not be used"), ns="httr")
f <- .prism_acis_cache_file(file.path(a[2],"cache","acis"),30.5,-84.5,"2020-01-01","2020-12-31")
mat <- .prism_fetch_acis_point(30.5,-84.5,"2020-01-01","2020-12-31",f)
stopifnot(.prism_complete_acis_matrix(mat,"2020-01-01","2020-12-31"))
process_weather_prism(data.frame(ID="P1",LAT=30.5,LONG=-84.5),2020,2020,file.path(a[2],"r-wth"),"ID","LAT","LONG",1,file.path(a[2],"r.log"),file.path(a[2],"cache"))
writeLines(jsonlite::toJSON(list(request=.prism_acis_request(30.5,-84.5,"2020-01-01","2020-12-31"),data=mat),auto_unbox=TRUE,digits=NA),f)
''')
    rdir = ROOT / "R"
    if (rdir / "cropmodeldata/R").is_dir():
        rdir = rdir / "cropmodeldata/R"
    subprocess.run([rscript, "--vanilla", str(script), str(rdir), str(tmp_path)],
                   capture_output=True, text=True, check=True, timeout=30)
    assert (tmp_path / "r-wth/P1.WTH").read_bytes() == (tmp_path / "wth/P1.WTH").read_bytes()
    cache = prism._acis_cache_file(str(tmp_path / "cache/acis"),30.5,-84.5,"2020-01-01","2020-12-31")
    assert prism._fetch_acis_point(30.5,-84.5,"2020-01-01","2020-12-31",cache) == records()
    assert prism.requests.post.call_count == 1

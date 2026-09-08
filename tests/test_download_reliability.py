from pathlib import Path
import pandas as pd
import geopandas as gpd
import pytest
from dssatutils import soil_polaris
from dssatutils.provider_retry import provider_retry, bounded_map, ProviderConnectivityError
from dssatutils.soil_validation import soil_file_issue


def test_mixed_polaris_preserves_good_point(tmp_path, monkeypatch):
    rows = pd.read_csv(Path(__file__).parent / 'fixtures/polaris_mixed_profiles.csv', dtype={'ID': str})
    monkeypatch.setattr(soil_polaris, '_fetch_polaris', lambda *a: rows)
    g = gpd.GeoDataFrame({'ID': ['00000001', '00000002']},
        geometry=gpd.points_from_xy([-84.6, -84.61], [30.55, 30.55]), crs=4326)
    soil_polaris.process_soils_polaris(g, str(tmp_path/'map.csv'), str(tmp_path), 'ID')
    assert soil_file_issue(tmp_path/'00000001.SOL') is None
    assert not (tmp_path/'00000002.SOL').exists()
    assert '00000002' in (tmp_path/'soil_processing_errors.log').read_text()


def test_network_recovery_is_bounded_and_preserves_errors():
    sleeps = []; attempts = []
    def fail():
        attempts.append(1)
        raise RuntimeError('Could not resolve hostname [example.org]')
    with pytest.raises(ProviderConnectivityError):
        provider_retry(fail, sleep=sleeps.append)
    assert len(attempts) == 3 and sleeps == [5, 10]
    calls = []
    def ordinary():
        calls.append(1)
        raise ValueError('invalid response schema')
    with pytest.raises(ValueError):
        provider_retry(ordinary, sleep=sleeps.append)
    assert len(calls) == 1


def test_stopped_network_does_not_consume_whole_queue():
    calls = []
    def fail(job):
        calls.append(job)
        raise ProviderConnectivityError('offline')
    with pytest.raises(ProviderConnectivityError):
        list(bounded_map(fail, range(10000), 1))
    assert calls == [0]


def test_missing_hydraulics_are_rejected(tmp_path, monkeypatch):
    test_mixed_polaris_preserves_good_point(tmp_path, monkeypatch)
    path = tmp_path/'00000001.SOL'
    lines = path.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith('@  SLB')) + 1
    row = lines[index]
    # SDUL occupies columns 20..24 in the DSSAT fixed-width layer table.
    import re
    header = lines[index-1]
    span = next(m.span() for m in re.finditer(r'\S+', header) if m.group() == 'SDUL')
    end = span[1]; start = end-5
    lines[index] = row[:start] + '  -99' + row[end:]
    path.write_text('\n'.join(lines)+'\n')
    assert 'hydraulic' in soil_file_issue(path)

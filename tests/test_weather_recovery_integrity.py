"""Offline regressions for annual completeness and real cache ownership."""
import os
import shutil
import subprocess
import sys
import threading
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'python'))
from dssatutils import weather_agera5 as ag
from dssatutils.weather_validation import is_wth_valid

FIXTURE = Path(__file__).parent / 'fixtures/agera5_2001_complete.csv'
AREA = [40.05, -90.05, 39.95, -89.95]


def test_requested_calendar_including_leap_and_explicit_cutoff(tmp_path):
    path = tmp_path / 'x.WTH'
    dates = pd.date_range('2000-01-01', '2000-12-31')
    rows = [f'{d:%Y%j} 12 25 15 2 10 60 3' for d in dates]
    def write(items): path.write_text('\n'.join(items) + '\n')
    write(rows)
    assert is_wth_valid(path, 2000, start_year=2000)
    for truncated in [rows[1:], rows[:-1], rows[182:183], rows[:59] + rows[60:]]:
        write(truncated)
        assert not is_wth_valid(path, 2000, start_year=2000)
    write(rows[:60])
    assert is_wth_valid(path, 2000, start_year=2000, end_date='2000-02-29')
    assert not is_wth_valid(path, 2000, start_year=2000, end_date='2000-03-01')


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -99])
def test_cache_rejects_single_missing_forcing(tmp_path, bad):
    x = pd.read_csv(FIXTURE)
    x['Solar_Radiation_Flux'] = x['Solar_Radiation_Flux'].astype(float)
    x.loc[5, 'Solar_Radiation_Flux'] = bad
    p = tmp_path / 'bad.csv'; x.to_csv(p, index=False)
    assert not ag._valid_agera5_timeseries_csv(str(p), 2001)


def test_cache_accepts_multicell_but_rejects_wrong_area(tmp_path):
    x = pd.read_csv(FIXTURE); other = x.copy(); other['longitude'] = -89.9
    p = tmp_path / 'tile.csv'; pd.concat([x, other]).to_csv(p, index=False)
    assert ag._valid_agera5_timeseries_csv(str(p), 2001)
    assert not ag._valid_agera5_timeseries_csv(str(p), 2001, area=AREA)
    pd.concat([x, other.iloc[:-1]]).to_csv(p, index=False)
    assert not ag._valid_agera5_timeseries_csv(str(p), 2001)


def test_download_checks_staging_after_nonpath_or_error_and_reuses(tmp_path, monkeypatch):
    calls = []
    class Client:
        def retrieve(self, dataset, req, target):
            calls.append(target)
            assert Path(target).parent != tmp_path
            shutil.copyfile(FIXTURE, target)
            raise RuntimeError('client error after complete transfer')
    monkeypatch.setattr(ag, '_make_cds_client', lambda _: Client())
    result = ag._download_agera5_timeseries(2001, AREA, str(tmp_path))
    assert result and ag._valid_agera5_timeseries_csv(result, 2001)
    assert ag._download_agera5_timeseries(2001, AREA, str(tmp_path)) == result
    assert len(calls) == 1
    assert not list(tmp_path.glob('.agera5-request-*'))


def test_cache_only_and_failed_transfer_preserve_original(tmp_path, monkeypatch):
    dest = Path(ag._agera5_timeseries_cache_path(str(tmp_path), 2001, AREA, 'csv'))
    dest.write_text('incomplete evidence\n')
    class Client:
        def retrieve(self, dataset, req, target):
            Path(target).write_text('bad response')
    monkeypatch.setattr(ag, '_make_cds_client', lambda _: Client())
    assert ag._download_agera5_timeseries(2001, AREA, str(tmp_path), cache_only=True) is None
    assert dest.read_text() == 'incomplete evidence\n'
    assert ag._download_agera5_timeseries(2001, AREA, str(tmp_path)) is None
    assert dest.read_text() == 'incomplete evidence\n'


def _has_r_filelock() -> bool:
    if not shutil.which("Rscript"):
        return False
    try:
        res = subprocess.run(
            ["Rscript", "--vanilla", "-e", "stopifnot(requireNamespace('filelock', quietly = TRUE))"],
            capture_output=True,
            timeout=10,
        )
        return res.returncode == 0
    except Exception:
        return False


def test_thread_and_r_process_contend_for_same_lock(tmp_path):
    path = str(tmp_path / 'cache.lock')
    owner = ag._agera5_acquire_lock(path)
    try:
        os.utime(path, (1, 1))  # age must not invalidate ownership
        attempts = []
        thread = threading.Thread(target=lambda: attempts.append(ag._agera5_acquire_lock(path, .05)))
        thread.start(); thread.join(2)
        assert attempts == [None]
        if _has_r_filelock():
            proc = subprocess.run(['Rscript', '--vanilla', '-e',
                'a<-commandArgs(TRUE); l<-filelock::lock(a[1],timeout=50); stopifnot(is.null(l))', path],
                capture_output=True, text=True, timeout=30)
            assert proc.returncode == 0, proc.stderr
    finally:
        ag._agera5_release_lock(owner)
    assert Path(path).exists()  # persistent inode, not a stale active lock
    second = ag._agera5_acquire_lock(path, .1)
    assert second is not None
    ag._agera5_release_lock(second)


def test_r_owner_blocks_python_and_crash_releases_lock(tmp_path):
    if not _has_r_filelock():
        pytest.skip("Rscript with filelock package is required for R cross-language lock contention test")
    path = str(tmp_path / 'cache.lock')
    proc = subprocess.Popen(['Rscript', '--vanilla', '-e',
        'a<-commandArgs(TRUE); l<-filelock::lock(a[1]); cat("READY\\n"); flush(stdout()); Sys.sleep(30)', path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert proc.stdout.readline().strip() == 'READY'
        assert ag._agera5_acquire_lock(path, .05) is None
    finally:
        proc.kill(); proc.wait(timeout=5)
    lock = ag._agera5_acquire_lock(path, .1)
    assert lock is not None
    ag._agera5_release_lock(lock)


def test_python_process_crash_releases_lock(tmp_path):
    path = str(tmp_path / "cache_py.lock")
    code = (
        "import sys, time; from dssatutils import weather_agera5 as ag; "
        f"lock = ag._agera5_acquire_lock({path!r}); sys.stdout.write('READY\\n'); sys.stdout.flush(); time.sleep(30)"
    )
    proc = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        assert proc.stdout.readline().strip() == "READY"
        assert ag._agera5_acquire_lock(path, 0.05) is None
    finally:
        proc.kill()
        proc.wait(timeout=5)
    lock = ag._agera5_acquire_lock(path, 0.1)
    assert lock is not None
    ag._agera5_release_lock(lock)


"""Unit tests for centralized executable discovery."""
from __future__ import annotations

import os
import sys
import pytest

import dssatutils
from dssatutils.discovery import find_rscript, find_dssat, find_mpi_runner


def test_discovery_exports_on_package():
    assert hasattr(dssatutils, "find_rscript")
    assert hasattr(dssatutils, "find_dssat")
    assert hasattr(dssatutils, "find_mpi_runner")
    assert callable(dssatutils.find_rscript)
    assert callable(dssatutils.find_dssat)
    assert callable(dssatutils.find_mpi_runner)


def test_find_rscript_explicit_nonexistent():
    assert find_rscript("non_existent_binary_xyz_12345") is None


def test_find_rscript_explicit_valid():
    res = find_rscript(sys.executable)
    assert res is not None
    assert os.path.isabs(res)
    assert os.path.isfile(res)


def test_find_rscript_env_var(monkeypatch, tmp_path):
    fake_exe = tmp_path / ("fake_rscript.exe" if sys.platform == "win32" else "fake_rscript")
    fake_exe.write_text("#!/bin/sh\necho 1\n")
    fake_exe.chmod(0o755)
    monkeypatch.setenv("RSCRIPT", str(fake_exe))
    res = find_rscript()
    assert res == str(fake_exe.resolve())


def test_find_dssat_env_var(monkeypatch, tmp_path):
    fake_exe = tmp_path / ("dscsm048.exe" if sys.platform == "win32" else "dscsm048")
    fake_exe.write_text("#!/bin/sh\necho DSSAT\n")
    fake_exe.chmod(0o755)
    monkeypatch.setenv("DSSAT_EXE", str(fake_exe))
    res = find_dssat()
    assert res == str(fake_exe.resolve())


def test_find_mpi_runner_env_var(monkeypatch, tmp_path):
    fake_exe = tmp_path / ("mpiexec.exe" if sys.platform == "win32" else "mpiexec")
    fake_exe.write_text("#!/bin/sh\necho MPI\n")
    fake_exe.chmod(0o755)
    monkeypatch.setenv("MPIEXEC", str(fake_exe))
    res = find_mpi_runner()
    assert res == str(fake_exe.resolve())


def test_find_rscript_system():
    # If Rscript is available on this system, find_rscript should return an absolute path
    import shutil
    has_system_r = shutil.which("Rscript")
    res = find_rscript()
    if has_system_r:
        assert res is not None
        assert os.path.isabs(res)
        assert os.path.exists(res)

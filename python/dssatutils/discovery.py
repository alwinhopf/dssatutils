"""Centralized executable discovery for DSSAT pipelines.

Resolves absolute paths for external tools (Rscript, DSSAT CSM, MPI runners)
to avoid bare executable invocations in subprocesses.
"""
from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path
from typing import Optional, Union

PathLike = Union[str, os.PathLike]

__all__ = ["find_rscript", "find_dssat", "find_mpi_runner"]


def _is_executable_file(path: str) -> bool:
    """Return True if path points to an existing regular file and is executable."""
    if not os.path.isfile(path):
        return False
    if sys.platform == "win32":
        return True
    return os.access(path, os.X_OK)


def _resolve_candidate(candidate: PathLike) -> Optional[str]:
    """Given a candidate string or path, resolve to an absolute executable path if valid."""
    str_cand = os.fspath(candidate).strip()
    if not str_cand:
        return None

    # Check if candidate has directory components or exists as an absolute/relative path
    has_dir = os.path.sep in str_cand or (os.path.altsep and os.path.altsep in str_cand)
    expanded = os.path.abspath(os.path.expanduser(str_cand))

    if has_dir or os.path.exists(expanded):
        if _is_executable_file(expanded):
            return expanded
        if sys.platform == "win32":
            for ext in (".exe", ".cmd", ".bat"):
                with_ext = expanded + ext
                if _is_executable_file(with_ext):
                    return with_ext

    # Search in PATH
    found = shutil.which(str_cand)
    if found and _is_executable_file(found):
        return os.path.abspath(found)

    if sys.platform == "win32" and not str_cand.lower().endswith((".exe", ".cmd", ".bat")):
        for ext in (".exe", ".cmd", ".bat"):
            found = shutil.which(str_cand + ext)
            if found and _is_executable_file(found):
                return os.path.abspath(found)

    return None


def find_rscript(explicit: Optional[PathLike] = None) -> Optional[str]:
    """Find the Rscript executable, returning an absolute path or None.

    Search order:
    1. `explicit` argument (if supplied).
    2. `RSCRIPT` environment variable.
    3. Standard `shutil.which` lookups ('Rscript', 'Rscript.exe').
    4. Common R installation directories (macOS Frameworks, Linux /usr/bin, Windows Program Files).
    """
    if explicit is not None and os.fspath(explicit).strip():
        return _resolve_candidate(explicit)

    env_val = os.environ.get("RSCRIPT", "").strip()
    if env_val:
        resolved = _resolve_candidate(env_val)
        if resolved:
            return resolved

    candidates = ["Rscript", "Rscript.exe"]
    for cand in candidates:
        resolved = _resolve_candidate(cand)
        if resolved:
            return resolved

    # Platform-specific fallback locations
    known_locations: list[str] = []
    if sys.platform == "darwin":
        known_locations.extend([
            "/usr/local/bin/Rscript",
            "/opt/homebrew/bin/Rscript",
            "/Library/Frameworks/R.framework/Resources/bin/Rscript",
        ])
    elif sys.platform == "win32":
        for pf in (os.environ.get("ProgramFiles", r"C:\Program Files"),
                   os.environ.get("ProgramFiles(x86)", r"C:\Program Files (x86)")):
            r_dir = Path(pf) / "R"
            if r_dir.is_dir():
                for sub in sorted(r_dir.glob("R-*"), reverse=True):
                    cand_path = sub / "bin" / "Rscript.exe"
                    known_locations.append(str(cand_path))
                    cand_x64 = sub / "bin" / "x64" / "Rscript.exe"
                    known_locations.append(str(cand_x64))
    else:  # Linux / Unix
        known_locations.extend([
            "/usr/bin/Rscript",
            "/usr/local/bin/Rscript",
        ])

    for loc in known_locations:
        if _is_executable_file(loc):
            return os.path.abspath(loc)

    return None


def find_dssat(explicit: Optional[PathLike] = None) -> Optional[str]:
    """Find the DSSAT CSM executable (e.g. dscsm048), returning an absolute path or None.

    Search order:
    1. `explicit` argument (if supplied).
    2. `DSSAT_EXE` environment variable.
    3. `DSSAT_HOME` or `DSSAT_DIR` environment variables combined with binary names.
    4. PATH search for standard versioned and unversioned names (`dscsm048`, `dscsm047`, `dscsm`).
    5. Common DSSAT install paths (/opt/dssat48, C:/DSSAT48, etc.).
    """
    if explicit is not None and os.fspath(explicit).strip():
        return _resolve_candidate(explicit)

    env_val = os.environ.get("DSSAT_EXE", "").strip()
    if env_val:
        resolved = _resolve_candidate(env_val)
        if resolved:
            return resolved

    binary_names = [
        "dscsm048", "dscsm048.exe",
        "dscsm047", "dscsm047.exe",
        "dscsm046", "dscsm046.exe",
        "dscsm045", "dscsm045.exe",
        "dscsm", "dscsm.exe",
    ]

    for env_dir_key in ("DSSAT_HOME", "DSSAT_DIR"):
        env_dir = os.environ.get(env_dir_key, "").strip()
        if env_dir and os.path.isdir(env_dir):
            for name in binary_names:
                candidate = os.path.join(env_dir, name)
                if _is_executable_file(candidate):
                    return os.path.abspath(candidate)

    for name in binary_names:
        resolved = _resolve_candidate(name)
        if resolved:
            return resolved

    # Common DSSAT installation paths
    search_dirs: list[str] = []
    if sys.platform == "win32":
        search_dirs.extend([r"C:\DSSAT48", r"C:\DSSAT47", r"C:\DSSAT46", r"C:\dssat48"])
    else:
        search_dirs.extend([
            "/opt/dssat48",
            "/usr/local/dssat48",
            os.path.expanduser("~/DSSAT48"),
            os.path.expanduser("~/dssat48"),
        ])

    for sdir in search_dirs:
        if os.path.isdir(sdir):
            for name in binary_names:
                candidate = os.path.join(sdir, name)
                if _is_executable_file(candidate):
                    return os.path.abspath(candidate)

    return None


def find_mpi_runner(explicit: Optional[PathLike] = None) -> Optional[str]:
    """Find an MPI runner executable (mpiexec, mpirun, srun), returning an absolute path or None.

    Search order:
    1. `explicit` argument (if supplied).
    2. `MPIEXEC` or `MPIRUN` environment variables.
    3. PATH search for standard MPI runners ('mpiexec', 'mpirun', 'srun').
    """
    if explicit is not None and os.fspath(explicit).strip():
        return _resolve_candidate(explicit)

    for env_key in ("MPIEXEC", "MPIRUN"):
        env_val = os.environ.get(env_key, "").strip()
        if env_val:
            resolved = _resolve_candidate(env_val)
            if resolved:
                return resolved

    mpi_names = [
        "mpiexec", "mpiexec.exe",
        "mpirun", "mpirun.exe",
        "srun",
    ]
    for name in mpi_names:
        resolved = _resolve_candidate(name)
        if resolved:
            return resolved

    return None

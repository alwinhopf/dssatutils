"""General DSSAT .WTH weather repair utilities.

Short gaps and invalid temperature pairs are fixed after weather retrieval, not
inside individual providers, so optional behavior applies to every weather
source.
"""

from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

DEFAULT_WEATHER_REPAIR_VARS = ("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")


def _append_log(path: str | Path | None, lines: Iterable[str]) -> None:
    if not path:
        return
    lines = list(lines)
    if not lines:
        return
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("a", encoding="utf-8") as fh:
        for line in lines:
            fh.write(str(line) + "\n")


def _timestamp() -> str:
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _find_header(lines: list[str]) -> int | None:
    for i, line in enumerate(lines):
        parts = line.strip().upper().split()
        if parts == ["@", "DATE", "SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND"]:
            return i
    return None


def _date_label(value) -> str:
    d = _date_from_code(value)
    if d is not None:
        return str(d)
    return str(value).strip()


def _date_from_code(value):
    txt = str(value).strip()
    if len(txt) < 4 or not txt.isdigit():
        return None
    try:
        year = int(txt[:-3])
        doy = int(txt[-3:])
        if year < 100:
            year = 1900 + year if year >= 80 else 2000 + year
        out = datetime(year, 1, 1).date() + timedelta(days=doy - 1)
        if out.year != year:
            return None
        return out
    except Exception:  # noqa: BLE001
        return None


def _code_from_date(value) -> str:
    return f"{value.year}{value.timetuple().tm_yday:03d}"


def _parse_wth(path: str | Path, log_file: str | Path | None = None, issue: str = "WEATHER_QA"):
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Weather file not found: {path}")
    lines = path.read_text(encoding="utf-8").splitlines()
    header_idx = _find_header(lines)
    pid = path.stem
    if header_idx is None or header_idx >= len(lines) - 1:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue={issue} "
            "status=skipped reason=no_DSSAT_weather_header"
        ])
        return lines, header_idx, pid, pd.DataFrame(), "skipped_no_header"

    data_lines = [line for line in lines[header_idx + 1:] if line.strip() and line.strip()[0].isdigit()]
    if not data_lines:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue={issue} "
            "status=skipped reason=no_weather_rows"
        ])
        return lines, header_idx, pid, pd.DataFrame(), "skipped_no_rows"

    rows = []
    for line in data_lines:
        parts = line.split()
        if len(parts) >= 8:
            rows.append(parts[:8])
    dat = pd.DataFrame(rows, columns=("DATE", *DEFAULT_WEATHER_REPAIR_VARS))
    for col in DEFAULT_WEATHER_REPAIR_VARS:
        dat[col] = pd.to_numeric(dat[col], errors="coerce")
        dat.loc[np.isclose(dat[col], -99.0), col] = np.nan
    dat["_DATE_OBJ"] = dat["DATE"].map(_date_from_code)
    return lines, header_idx, pid, dat, "ok"


def _write_wth(path: str | Path, lines: list[str], header_idx: int, dat: pd.DataFrame) -> None:
    write_dat = dat.copy()
    if "_DATE_OBJ" in write_dat.columns:
        write_dat = write_dat.drop(columns=["_DATE_OBJ"])
    for col in DEFAULT_WEATHER_REPAIR_VARS:
        values = write_dat[col].to_numpy(dtype=float).copy()
        values[~np.isfinite(values)] = -99.0
        values[(values >= 9999.95) | (values <= -999.95)] = -99.0
        write_dat[col] = values
    formatted = []
    for _, row in write_dat.iterrows():
        line = (
            f"{str(row['DATE']):>7s}"
            f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
            f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}"
            f"{row['WIND']:6.1f}"
        )
        formatted.append(line.replace("-99.0", "  -99"))
    Path(path).write_text("\n".join(lines[:header_idx + 1] + formatted) + "\n", encoding="utf-8")


def _missing_runs(mask: np.ndarray) -> list[tuple[int, int, int]]:
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return []
    runs: list[tuple[int, int, int]] = []
    start = int(idx[0])
    prev = int(idx[0])
    for value in idx[1:]:
        value = int(value)
        if value == prev + 1:
            prev = value
        else:
            runs.append((start, prev, prev - start + 1))
            start = prev = value
    runs.append((start, prev, prev - start + 1))
    return runs


def repair_weather_file_missing_values(
    wth_file: str | Path,
    max_gap_days: int = 3,
    window_days: int = 2,
    variables: Iterable[str] = DEFAULT_WEATHER_REPAIR_VARS,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Repair short missing-value gaps in one DSSAT weather file.

    Missing values encoded as ``-99``/``-99.0`` are repaired only when the
    contiguous gap is no longer than ``max_gap_days`` and the full two-day window
    before and after the gap is present. Each corrected value is filled with the
    mean of those surrounding days for the same weather variable.
    """

    path = Path(wth_file)
    if not path.exists():
        raise FileNotFoundError(f"Weather file not found: {path}")
    max_gap_days = int(max_gap_days)
    window_days = int(window_days)
    if max_gap_days < 1:
        raise ValueError("max_gap_days must be >= 1")
    if window_days < 1:
        raise ValueError("window_days must be >= 1")

    lines = path.read_text(encoding="utf-8").splitlines()
    header_idx = _find_header(lines)
    pid = path.stem
    if header_idx is None or header_idx >= len(lines) - 1:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} status=skipped reason=no_DSSAT_weather_header"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "variable": None,
            "repaired_count": 0, "unrepaired_count": 0, "repaired_runs": 0,
            "status": "skipped_no_header",
        }])

    data_lines = [line for line in lines[header_idx + 1:] if line.strip() and line.strip()[0].isdigit()]
    if not data_lines:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} status=skipped reason=no_weather_rows"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "variable": None,
            "repaired_count": 0, "unrepaired_count": 0, "repaired_runs": 0,
            "status": "skipped_no_rows",
        }])

    rows = []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 8:
            continue
        rows.append(parts[:8])
    dat = pd.DataFrame(rows, columns=("DATE", *DEFAULT_WEATHER_REPAIR_VARS))
    for col in DEFAULT_WEATHER_REPAIR_VARS:
        dat[col] = pd.to_numeric(dat[col], errors="coerce")
        dat.loc[np.isclose(dat[col], -99.0), col] = np.nan

    wanted = [v.upper() for v in variables]
    wanted = [v for v in wanted if v in DEFAULT_WEATHER_REPAIR_VARS]
    if not wanted:
        wanted = list(DEFAULT_WEATHER_REPAIR_VARS)

    original = dat.copy(deep=True)
    log_lines: list[str] = []
    summary: list[dict] = []

    for var in wanted:
        before = original[var].to_numpy(dtype=float)
        missing = ~np.isfinite(before)
        repaired_count = 0
        unrepaired_count = 0
        repaired_runs = 0
        for start, end, length in _missing_runs(missing):
            if length <= max_gap_days:
                neighbor_idx = list(range(start - window_days, start)) + list(range(end + 1, end + 1 + window_days))
                in_bounds = all(0 <= i < len(original) for i in neighbor_idx)
                neighbor_vals = before[neighbor_idx] if in_bounds else np.array([])
                usable = in_bounds and len(neighbor_vals) == 2 * window_days and np.all(np.isfinite(neighbor_vals))
                if usable:
                    fill_value = float(np.mean(neighbor_vals))
                    dat.loc[start:end, var] = fill_value
                    repaired_count += length
                    repaired_runs += 1
                    log_lines.append(
                        f"{_timestamp()} file={path.name} id={pid} variable={var} status=repaired "
                        f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                        f"gap_days={length} fill_value={fill_value:.4f} method=mean_{window_days}_days_before_after "
                        f"neighbor_dates={_date_label(dat.loc[neighbor_idx[0], 'DATE'])}.."
                        f"{_date_label(dat.loc[neighbor_idx[window_days - 1], 'DATE'])};"
                        f"{_date_label(dat.loc[neighbor_idx[window_days], 'DATE'])}.."
                        f"{_date_label(dat.loc[neighbor_idx[-1], 'DATE'])}"
                    )
                else:
                    unrepaired_count += length
                    log_lines.append(
                        f"{_timestamp()} file={path.name} id={pid} variable={var} status=unrepaired "
                        f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                        f"gap_days={length} reason=insufficient_{window_days}_day_neighbors"
                    )
            else:
                unrepaired_count += length
                log_lines.append(
                    f"{_timestamp()} file={path.name} id={pid} variable={var} status=unrepaired "
                    f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                    f"gap_days={length} reason=gap_exceeds_max_{max_gap_days}_days"
                )

        summary.append({
            "file": str(path),
            "id": pid,
            "variable": var,
            "repaired_count": repaired_count,
            "unrepaired_count": unrepaired_count,
            "repaired_runs": repaired_runs,
            "status": "repaired" if repaired_count else ("unrepaired_missing" if unrepaired_count else "unchanged"),
        })

    _append_log(log_file, log_lines)

    if not dry_run and any(row["repaired_count"] > 0 for row in summary):
        write_dat = dat.copy()
        for col in DEFAULT_WEATHER_REPAIR_VARS:
            values = write_dat[col].to_numpy(dtype=float).copy()
            values[~np.isfinite(values)] = -99.0
            values[(values >= 9999.95) | (values <= -999.95)] = -99.0
            write_dat[col] = values
        formatted = []
        for _, row in write_dat.iterrows():
            line = (
                f"{str(row['DATE']):>7s}"
                f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
                f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}"
                f"{row['WIND']:6.1f}"
            )
            formatted.append(line.replace("-99.0", "  -99"))
        path.write_text("\n".join(lines[:header_idx + 1] + formatted) + "\n", encoding="utf-8")

    return pd.DataFrame(summary)


def repair_weather_missing_values(
    weather_dir: str | Path,
    ids: Iterable[str] | None = None,
    max_gap_days: int = 3,
    window_days: int = 2,
    variables: Iterable[str] = DEFAULT_WEATHER_REPAIR_VARS,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Repair short missing-value gaps in DSSAT weather files."""

    weather_dir = Path(weather_dir)
    if not weather_dir.exists():
        raise FileNotFoundError(f"Weather directory not found: {weather_dir}")
    if log_file is None:
        log_file = weather_dir / "weather_repair.log"
    files = sorted(weather_dir.glob("*.WTH"))
    if ids:
        wanted = {f"{str(i)}.WTH" for i in ids}
        files = [f for f in files if f.name in wanted]
    if not files:
        _append_log(log_file, [
            f"{_timestamp()} weather_dir={weather_dir} status=skipped reason=no_weather_files"
        ])
        return pd.DataFrame(columns=[
            "file", "id", "variable", "repaired_count", "unrepaired_count", "repaired_runs", "status"
        ])

    _append_log(log_file, [
        "",
        f"{_timestamp()} weather_dir={weather_dir} status=started files={len(files)} "
        f"max_gap_days={int(max_gap_days)} window_days={int(window_days)} "
        f"variables={','.join(variables)} dry_run={dry_run}",
    ])

    parts = [
        repair_weather_file_missing_values(
            f,
            max_gap_days=max_gap_days,
            window_days=window_days,
            variables=variables,
            log_file=log_file,
            dry_run=dry_run,
        )
        for f in files
    ]
    summary = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    repaired = int(summary["repaired_count"].sum()) if "repaired_count" in summary else 0
    unrepaired = int(summary["unrepaired_count"].sum()) if "unrepaired_count" in summary else 0
    _append_log(log_file, [
        f"{_timestamp()} weather_dir={weather_dir} status=finished "
        f"repaired_values={repaired} unrepaired_values={unrepaired} log_file={log_file}"
    ])
    return summary


def repair_weather_file_temperature_inversions(
    wth_file: str | Path,
    max_gap_days: int = 3,
    window_days: int = 2,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Repair short Tmax/Tmin inversion runs in one DSSAT weather file.

    Rows where both temperatures are present but ``TMIN > TMAX`` are repaired
    only when the contiguous inversion run is no longer than ``max_gap_days`` and
    the full ``window_days`` before and after the run has valid, non-inverted
    temperatures. ``TMAX`` and ``TMIN`` are replaced independently using the mean
    of their respective neighboring values.
    """

    path = Path(wth_file)
    if not path.exists():
        raise FileNotFoundError(f"Weather file not found: {path}")
    max_gap_days = int(max_gap_days)
    window_days = int(window_days)
    if max_gap_days < 1:
        raise ValueError("max_gap_days must be >= 1")
    if window_days < 1:
        raise ValueError("window_days must be >= 1")

    lines = path.read_text(encoding="utf-8").splitlines()
    header_idx = _find_header(lines)
    pid = path.stem
    if header_idx is None or header_idx >= len(lines) - 1:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue=TMIN_GT_TMAX "
            "status=skipped reason=no_DSSAT_weather_header"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "issue": "TMIN_GT_TMAX",
            "repaired_count": 0, "unrepaired_count": 0, "repaired_runs": 0,
            "status": "skipped_no_header",
        }])

    data_lines = [line for line in lines[header_idx + 1:] if line.strip() and line.strip()[0].isdigit()]
    if not data_lines:
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue=TMIN_GT_TMAX "
            "status=skipped reason=no_weather_rows"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "issue": "TMIN_GT_TMAX",
            "repaired_count": 0, "unrepaired_count": 0, "repaired_runs": 0,
            "status": "skipped_no_rows",
        }])

    rows = []
    for line in data_lines:
        parts = line.split()
        if len(parts) < 8:
            continue
        rows.append(parts[:8])
    dat = pd.DataFrame(rows, columns=("DATE", *DEFAULT_WEATHER_REPAIR_VARS))
    for col in DEFAULT_WEATHER_REPAIR_VARS:
        dat[col] = pd.to_numeric(dat[col], errors="coerce")
        dat.loc[np.isclose(dat[col], -99.0), col] = np.nan

    original = dat.copy(deep=True)
    tmax = original["TMAX"].to_numpy(dtype=float)
    tmin = original["TMIN"].to_numpy(dtype=float)
    inversion = np.isfinite(tmax) & np.isfinite(tmin) & (tmin > tmax)

    log_lines: list[str] = []
    repaired_count = 0
    unrepaired_count = 0
    repaired_runs = 0

    for start, end, length in _missing_runs(inversion):
        if length <= max_gap_days:
            neighbor_idx = list(range(start - window_days, start)) + list(range(end + 1, end + 1 + window_days))
            in_bounds = all(0 <= i < len(original) for i in neighbor_idx)
            neighbor_tmax = tmax[neighbor_idx] if in_bounds else np.array([])
            neighbor_tmin = tmin[neighbor_idx] if in_bounds else np.array([])
            usable = (
                in_bounds
                and len(neighbor_tmax) == 2 * window_days
                and np.all(np.isfinite(neighbor_tmax))
                and np.all(np.isfinite(neighbor_tmin))
                and np.all(neighbor_tmin <= neighbor_tmax)
            )
            if usable:
                fill_tmax = float(np.mean(neighbor_tmax))
                fill_tmin = float(np.mean(neighbor_tmin))
                dat.loc[start:end, "TMAX"] = fill_tmax
                dat.loc[start:end, "TMIN"] = fill_tmin
                repaired_count += length
                repaired_runs += 1
                log_lines.append(
                    f"{_timestamp()} file={path.name} id={pid} issue=TMIN_GT_TMAX status=repaired "
                    f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                    f"gap_days={length} fill_TMAX={fill_tmax:.4f} fill_TMIN={fill_tmin:.4f} "
                    f"method=mean_{window_days}_days_before_after "
                    f"neighbor_dates={_date_label(dat.loc[neighbor_idx[0], 'DATE'])}.."
                    f"{_date_label(dat.loc[neighbor_idx[window_days - 1], 'DATE'])};"
                    f"{_date_label(dat.loc[neighbor_idx[window_days], 'DATE'])}.."
                    f"{_date_label(dat.loc[neighbor_idx[-1], 'DATE'])}"
                )
            else:
                unrepaired_count += length
                log_lines.append(
                    f"{_timestamp()} file={path.name} id={pid} issue=TMIN_GT_TMAX status=unrepaired "
                    f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                    f"gap_days={length} reason=insufficient_{window_days}_day_valid_temperature_neighbors"
                )
        else:
            unrepaired_count += length
            log_lines.append(
                f"{_timestamp()} file={path.name} id={pid} issue=TMIN_GT_TMAX status=unrepaired "
                f"dates={_date_label(dat.loc[start, 'DATE'])}..{_date_label(dat.loc[end, 'DATE'])} "
                f"gap_days={length} reason=gap_exceeds_max_{max_gap_days}_days"
            )

    _append_log(log_file, log_lines)

    if not dry_run and repaired_count > 0:
        write_dat = dat.copy()
        for col in DEFAULT_WEATHER_REPAIR_VARS:
            values = write_dat[col].to_numpy(dtype=float).copy()
            values[~np.isfinite(values)] = -99.0
            values[(values >= 9999.95) | (values <= -999.95)] = -99.0
            write_dat[col] = values
        formatted = []
        for _, row in write_dat.iterrows():
            line = (
                f"{str(row['DATE']):>7s}"
                f"{row['SRAD']:6.1f}{row['TMAX']:6.1f}{row['TMIN']:6.1f}"
                f"{row['RAIN']:6.1f}{row['TDEW']:6.1f}{row['RH2M']:6.1f}"
                f"{row['WIND']:6.1f}"
            )
            formatted.append(line.replace("-99.0", "  -99"))
        path.write_text("\n".join(lines[:header_idx + 1] + formatted) + "\n", encoding="utf-8")

    return pd.DataFrame([{
        "file": str(path),
        "id": pid,
        "issue": "TMIN_GT_TMAX",
        "repaired_count": repaired_count,
        "unrepaired_count": unrepaired_count,
        "repaired_runs": repaired_runs,
        "status": "repaired" if repaired_count else (
            "unrepaired_temperature_inversion" if unrepaired_count else "unchanged"
        ),
    }])


def repair_weather_temperature_inversions(
    weather_dir: str | Path,
    ids: Iterable[str] | None = None,
    max_gap_days: int = 3,
    window_days: int = 2,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Repair short Tmax/Tmin inversion runs in DSSAT weather files."""

    weather_dir = Path(weather_dir)
    if not weather_dir.exists():
        raise FileNotFoundError(f"Weather directory not found: {weather_dir}")
    if log_file is None:
        log_file = weather_dir / "weather_repair.log"
    files = sorted(weather_dir.glob("*.WTH"))
    if ids:
        wanted = {f"{str(i)}.WTH" for i in ids}
        files = [f for f in files if f.name in wanted]
    if not files:
        _append_log(log_file, [
            f"{_timestamp()} weather_dir={weather_dir} issue=TMIN_GT_TMAX "
            "status=skipped reason=no_weather_files"
        ])
        return pd.DataFrame(columns=[
            "file", "id", "issue", "repaired_count", "unrepaired_count", "repaired_runs", "status"
        ])

    _append_log(log_file, [
        "",
        f"{_timestamp()} weather_dir={weather_dir} issue=TMIN_GT_TMAX status=started "
        f"files={len(files)} max_gap_days={int(max_gap_days)} "
        f"window_days={int(window_days)} dry_run={dry_run}",
    ])

    parts = [
        repair_weather_file_temperature_inversions(
            f,
            max_gap_days=max_gap_days,
            window_days=window_days,
            log_file=log_file,
            dry_run=dry_run,
        )
        for f in files
    ]
    summary = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    repaired = int(summary["repaired_count"].sum()) if "repaired_count" in summary else 0
    unrepaired = int(summary["unrepaired_count"].sum()) if "unrepaired_count" in summary else 0
    _append_log(log_file, [
        f"{_timestamp()} weather_dir={weather_dir} issue=TMIN_GT_TMAX status=finished "
        f"repaired_values={repaired} unrepaired_values={unrepaired} log_file={log_file}"
    ])
    return summary


def repair_weather_file_date_gaps(
    wth_file: str | Path,
    max_gap_days: int = 3,
    window_days: int = 2,
    variables: Iterable[str] = DEFAULT_WEATHER_REPAIR_VARS,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Insert short missing DATE rows in one DSSAT weather file.

    Missing calendar rows are repaired only when the contiguous missing-date run
    is no longer than ``max_gap_days`` and all ``window_days`` before and after
    the run exist with finite values. Inserted rows receive the mean of the
    neighboring values for each requested variable.
    """

    path = Path(wth_file)
    max_gap_days = int(max_gap_days)
    window_days = int(window_days)
    if max_gap_days < 1:
        raise ValueError("max_gap_days must be >= 1")
    if window_days < 1:
        raise ValueError("window_days must be >= 1")
    lines, header_idx, pid, dat, status = _parse_wth(path, log_file, issue="DATE_GAP")
    if status != "ok":
        return pd.DataFrame([{
            "file": str(path), "id": pid, "issue": "DATE_GAP",
            "repaired_count": 0, "unrepaired_count": 0, "repaired_runs": 0,
            "status": status,
        }])
    if dat["_DATE_OBJ"].isna().any():
        bad = int(dat["_DATE_OBJ"].isna().sum())
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue=DATE_GAP status=unrepaired "
            f"reason=unparseable_date_codes count={bad}"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "issue": "DATE_GAP",
            "repaired_count": 0, "unrepaired_count": bad, "repaired_runs": 0,
            "status": "unrepaired_unparseable_dates",
        }])
    dup = dat["_DATE_OBJ"].duplicated(keep=False)
    if dup.any():
        bad = int(dup.sum())
        _append_log(log_file, [
            f"{_timestamp()} file={path.name} id={pid} issue=DATE_GAP status=unrepaired "
            f"reason=duplicate_dates count={bad}"
        ])
        return pd.DataFrame([{
            "file": str(path), "id": pid, "issue": "DATE_GAP",
            "repaired_count": 0, "unrepaired_count": bad, "repaired_runs": 0,
            "status": "unrepaired_duplicate_dates",
        }])

    wanted = [v.upper() for v in variables]
    wanted = [v for v in wanted if v in DEFAULT_WEATHER_REPAIR_VARS]
    if not wanted:
        wanted = list(DEFAULT_WEATHER_REPAIR_VARS)

    original_order = dat["_DATE_OBJ"].tolist()
    dat = dat.sort_values("_DATE_OBJ").reset_index(drop=True)
    sorted_rows = original_order != dat["_DATE_OBJ"].tolist()
    present = set(dat["_DATE_OBJ"].tolist())
    start_date = min(present)
    end_date = max(present)
    expected = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]
    missing = [d for d in expected if d not in present]

    runs = []
    if missing:
        run_start = run_prev = missing[0]
        for d in missing[1:]:
            if d == run_prev + timedelta(days=1):
                run_prev = d
            else:
                runs.append((run_start, run_prev, (run_prev - run_start).days + 1))
                run_start = run_prev = d
        runs.append((run_start, run_prev, (run_prev - run_start).days + 1))

    by_date = {d: i for i, d in enumerate(dat["_DATE_OBJ"].tolist())}
    new_rows = []
    log_lines: list[str] = []
    repaired_count = 0
    unrepaired_count = 0
    repaired_runs = 0

    for run_start, run_end, length in runs:
        if length <= max_gap_days:
            neighbor_dates = (
                [run_start - timedelta(days=i) for i in range(window_days, 0, -1)]
                + [run_end + timedelta(days=i) for i in range(1, window_days + 1)]
            )
            in_bounds = all(d in by_date for d in neighbor_dates)
            if in_bounds:
                neighbor_idx = [by_date[d] for d in neighbor_dates]
                usable = all(np.all(np.isfinite(dat.loc[neighbor_idx, var].to_numpy(dtype=float))) for var in wanted)
            else:
                neighbor_idx = []
                usable = False
            if usable:
                fill_values = {var: float(dat.loc[neighbor_idx, var].mean()) for var in wanted}
                for d in [run_start + timedelta(days=i) for i in range(length)]:
                    row = {"DATE": _code_from_date(d), "_DATE_OBJ": d}
                    for var in DEFAULT_WEATHER_REPAIR_VARS:
                        row[var] = fill_values.get(var, np.nan)
                    new_rows.append(row)
                repaired_count += length
                repaired_runs += 1
                fill_txt = ",".join(f"{k}={v:.4f}" for k, v in fill_values.items())
                log_lines.append(
                    f"{_timestamp()} file={path.name} id={pid} issue=DATE_GAP status=repaired "
                    f"dates={run_start}..{run_end} gap_days={length} fill_values={fill_txt} "
                    f"method=mean_{window_days}_days_before_after "
                    f"neighbor_dates={neighbor_dates[0]}..{neighbor_dates[window_days - 1]};"
                    f"{neighbor_dates[window_days]}..{neighbor_dates[-1]}"
                )
            else:
                unrepaired_count += length
                log_lines.append(
                    f"{_timestamp()} file={path.name} id={pid} issue=DATE_GAP status=unrepaired "
                    f"dates={run_start}..{run_end} gap_days={length} "
                    f"reason=insufficient_{window_days}_day_neighbors"
                )
        else:
            unrepaired_count += length
            log_lines.append(
                f"{_timestamp()} file={path.name} id={pid} issue=DATE_GAP status=unrepaired "
                f"dates={run_start}..{run_end} gap_days={length} reason=gap_exceeds_max_{max_gap_days}_days"
            )

    _append_log(log_file, log_lines)
    if not dry_run and (new_rows or sorted_rows):
        if new_rows:
            dat = pd.concat([dat, pd.DataFrame(new_rows)], ignore_index=True)
        dat = dat.sort_values("_DATE_OBJ").reset_index(drop=True)
        _write_wth(path, lines, int(header_idx), dat)

    status_out = "repaired" if repaired_count else ("sorted" if sorted_rows else ("unrepaired_date_gap" if unrepaired_count else "unchanged"))
    return pd.DataFrame([{
        "file": str(path), "id": pid, "issue": "DATE_GAP",
        "repaired_count": repaired_count, "unrepaired_count": unrepaired_count,
        "repaired_runs": repaired_runs, "status": status_out,
    }])


def repair_weather_date_gaps(
    weather_dir: str | Path,
    ids: Iterable[str] | None = None,
    max_gap_days: int = 3,
    window_days: int = 2,
    variables: Iterable[str] = DEFAULT_WEATHER_REPAIR_VARS,
    log_file: str | Path | None = None,
    dry_run: bool = False,
) -> pd.DataFrame:
    """Insert short missing DATE rows in DSSAT weather files."""

    weather_dir = Path(weather_dir)
    if not weather_dir.exists():
        raise FileNotFoundError(f"Weather directory not found: {weather_dir}")
    if log_file is None:
        log_file = weather_dir / "weather_repair.log"
    files = sorted(weather_dir.glob("*.WTH"))
    if ids:
        wanted = {f"{str(i)}.WTH" for i in ids}
        files = [f for f in files if f.name in wanted]
    if not files:
        _append_log(log_file, [
            f"{_timestamp()} weather_dir={weather_dir} issue=DATE_GAP status=skipped reason=no_weather_files"
        ])
        return pd.DataFrame(columns=[
            "file", "id", "issue", "repaired_count", "unrepaired_count", "repaired_runs", "status"
        ])

    _append_log(log_file, [
        "",
        f"{_timestamp()} weather_dir={weather_dir} issue=DATE_GAP status=started files={len(files)} "
        f"max_gap_days={int(max_gap_days)} window_days={int(window_days)} "
        f"variables={','.join(variables)} dry_run={dry_run}",
    ])
    parts = [
        repair_weather_file_date_gaps(
            f,
            max_gap_days=max_gap_days,
            window_days=window_days,
            variables=variables,
            log_file=log_file,
            dry_run=dry_run,
        )
        for f in files
    ]
    summary = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    repaired = int(summary["repaired_count"].sum()) if "repaired_count" in summary else 0
    unrepaired = int(summary["unrepaired_count"].sum()) if "unrepaired_count" in summary else 0
    _append_log(log_file, [
        f"{_timestamp()} weather_dir={weather_dir} issue=DATE_GAP status=finished "
        f"repaired_values={repaired} unrepaired_values={unrepaired} log_file={log_file}"
    ])
    return summary


def _add_issue(rows, path, pid, issue, severity, count, first_date="", last_date="", details=""):
    if int(count) <= 0:
        return
    rows.append({
        "file": str(path),
        "id": pid,
        "issue": issue,
        "severity": severity,
        "count": int(count),
        "first_date": str(first_date) if first_date else "",
        "last_date": str(last_date) if last_date else "",
        "details": details,
    })


def _flatline_runs(values: np.ndarray, dates: list, min_days: int):
    out = []
    start = 0
    while start < len(values):
        if not np.isfinite(values[start]):
            start += 1
            continue
        end = start
        while end + 1 < len(values) and np.isfinite(values[end + 1]) and values[end + 1] == values[start]:
            end += 1
        if end - start + 1 >= min_days:
            out.append((start, end, end - start + 1, values[start], dates[start], dates[end]))
        start = end + 1
    return out


def audit_weather_file_quality(
    wth_file: str | Path,
    flatline_days: int = 10,
    log_file: str | Path | None = None,
) -> pd.DataFrame:
    """Flag suspicious DSSAT weather rows without modifying the file."""

    path = Path(wth_file)
    lines, header_idx, pid, dat, status = _parse_wth(path, log_file, issue="WEATHER_QA")
    rows: list[dict] = []
    if status != "ok":
        _add_issue(rows, path, pid, status, "error", 1, details=status)
        return pd.DataFrame(rows)

    date_obj = dat["_DATE_OBJ"].tolist()
    bad_dates = [d is None or pd.isna(d) for d in date_obj]
    _add_issue(rows, path, pid, "unparseable_date_codes", "error", sum(bad_dates))
    if any(bad_dates):
        return pd.DataFrame(rows)

    dup = dat["_DATE_OBJ"].duplicated(keep=False)
    _add_issue(rows, path, pid, "duplicate_dates", "error", int(dup.sum()))
    out_of_order = any(date_obj[i] > date_obj[i + 1] for i in range(len(date_obj) - 1))
    _add_issue(rows, path, pid, "out_of_order_dates", "warning", int(out_of_order))

    unique_dates = sorted(set(date_obj))
    if unique_dates:
        expected = [unique_dates[0] + timedelta(days=i) for i in range((unique_dates[-1] - unique_dates[0]).days + 1)]
        missing_dates = [d for d in expected if d not in set(unique_dates)]
        _add_issue(
            rows, path, pid, "missing_date_rows", "error", len(missing_dates),
            first_date=missing_dates[0] if missing_dates else "",
            last_date=missing_dates[-1] if missing_dates else "",
        )

    for var in DEFAULT_WEATHER_REPAIR_VARS:
        missing = ~np.isfinite(dat[var].to_numpy(dtype=float))
        idx = np.flatnonzero(missing)
        _add_issue(
            rows, path, pid, f"{var}_missing_values", "warning", len(idx),
            first_date=_date_label(dat.loc[int(idx[0]), "DATE"]) if len(idx) else "",
            last_date=_date_label(dat.loc[int(idx[-1]), "DATE"]) if len(idx) else "",
        )

    tmax = dat["TMAX"].to_numpy(dtype=float)
    tmin = dat["TMIN"].to_numpy(dtype=float)
    rain = dat["RAIN"].to_numpy(dtype=float)
    srad = dat["SRAD"].to_numpy(dtype=float)
    rh = dat["RH2M"].to_numpy(dtype=float)
    wind = dat["WIND"].to_numpy(dtype=float)
    tdew = dat["TDEW"].to_numpy(dtype=float)

    checks = [
        ("tmin_gt_tmax", "error", np.isfinite(tmax) & np.isfinite(tmin) & (tmin > tmax), ""),
        ("tmax_out_of_range", "warning", np.isfinite(tmax) & ((tmax < -60) | (tmax > 60)), "bounds=-60..60C"),
        ("tmin_out_of_range", "warning", np.isfinite(tmin) & ((tmin < -70) | (tmin > 50)), "bounds=-70..50C"),
        ("diurnal_range_extreme", "warning", np.isfinite(tmax) & np.isfinite(tmin) & ((tmax - tmin) > 45), "TMAX-TMIN>45C"),
        ("rain_negative", "error", np.isfinite(rain) & (rain < 0), ""),
        ("rain_extreme", "warning", np.isfinite(rain) & (rain > 500), "RAIN>500mm"),
        ("srad_out_of_range", "warning", np.isfinite(srad) & ((srad < 0) | (srad > 40)), "bounds=0..40MJ/m2/day"),
        ("rh2m_out_of_range", "warning", np.isfinite(rh) & ((rh < 0) | (rh > 100)), "bounds=0..100%"),
        ("wind_out_of_range", "warning", np.isfinite(wind) & ((wind < 0) | (wind > 75)), "bounds=0..75m/s"),
        ("tdew_gt_tmax", "warning", np.isfinite(tdew) & np.isfinite(tmax) & (tdew > tmax), ""),
    ]
    for issue, severity, mask, details in checks:
        idx = np.flatnonzero(mask)
        _add_issue(
            rows, path, pid, issue, severity, len(idx),
            first_date=_date_label(dat.loc[int(idx[0]), "DATE"]) if len(idx) else "",
            last_date=_date_label(dat.loc[int(idx[-1]), "DATE"]) if len(idx) else "",
            details=details,
        )

    flat_vars = ("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
    for var in flat_vars:
        runs = _flatline_runs(dat[var].to_numpy(dtype=float), dat["_DATE_OBJ"].tolist(), int(flatline_days))
        if runs:
            total = sum(r[2] for r in runs)
            first = runs[0][4]
            last = runs[-1][5]
            values = ",".join(f"{r[3]:.2f}x{r[2]}d" for r in runs[:5])
            _add_issue(rows, path, pid, f"{var}_flatline", "info", total, first, last, f"min_days={flatline_days}; examples={values}")

    return pd.DataFrame(rows, columns=["file", "id", "issue", "severity", "count", "first_date", "last_date", "details"])


def audit_weather_quality(
    weather_dir: str | Path,
    ids: Iterable[str] | None = None,
    audit_csv: str | Path | None = None,
    flatline_days: int = 10,
    log_file: str | Path | None = None,
) -> pd.DataFrame:
    """Write a flag-only QA audit for DSSAT weather files."""

    weather_dir = Path(weather_dir)
    if not weather_dir.exists():
        raise FileNotFoundError(f"Weather directory not found: {weather_dir}")
    if log_file is None:
        log_file = weather_dir / "weather_repair.log"
    if audit_csv is None:
        audit_csv = weather_dir / "weather_quality_audit.csv"
    files = sorted(weather_dir.glob("*.WTH"))
    if ids:
        wanted = {f"{str(i)}.WTH" for i in ids}
        files = [f for f in files if f.name in wanted]
    parts = [audit_weather_file_quality(f, flatline_days=flatline_days, log_file=log_file) for f in files]
    summary = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(
        columns=["file", "id", "issue", "severity", "count", "first_date", "last_date", "details"]
    )
    Path(audit_csv).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(audit_csv, index=False)
    _append_log(log_file, [
        f"{_timestamp()} weather_dir={weather_dir} issue=WEATHER_QA status=finished "
        f"files={len(files)} findings={len(summary)} audit_csv={audit_csv}"
    ])
    return summary

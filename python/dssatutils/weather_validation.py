"""Structural validation for DSSAT fixed-width weather files."""
from __future__ import annotations

from datetime import date, timedelta
from pathlib import Path
import math
import re


def _parse_wth_data_line(line: str):
    fields = None
    if len(line) >= 49:
        fixed = [line[:7].strip()] + [line[start:start + 6].strip()
                                      for start in range(7, 49, 6)]
        if re.fullmatch(r"\d{5,7}", fixed[0]) and all(fixed[1:]):
            fields = fixed
    if fields is None:
        fields = line.split()
    if len(fields) != 8 or not re.fullmatch(r"\d{5,7}", fields[0]):
        return None
    try:
        values = [float(value) for value in fields[1:]]
    except ValueError:
        return None
    if not all(math.isfinite(value) for value in values):
        return None
    return fields[0], values


def _wth_code_to_date(code: str):
    if len(code) == 5:
        yy, doy = int(code[:2]), int(code[2:])
        year = 2000 + yy if yy < 80 else 1900 + yy
    elif len(code) == 7:
        year, doy = int(code[:4]), int(code[4:])
    else:
        return None
    try:
        value = date(year, 1, 1) + timedelta(days=doy - 1)
    except (OverflowError, ValueError):
        return None
    return value if doy >= 1 and value.year == year else None


def is_wth_valid(path: str | Path, end_year: int | None = None,
                 required_columns: tuple[str, ...] | list[str] | None = None) -> bool:
    """Return whether *path* is a complete, parseable DSSAT weather file.

    ``required_columns`` optionally rejects DSSAT ``-99`` missing markers in
    forcing variables required by the caller's model configuration.
    """
    target = Path(path)
    if not target.is_file() or target.stat().st_size <= 0:
        return False
    try:
        lines = target.read_text(encoding="utf-8").splitlines()
        parsed = [_parse_wth_data_line(line) for line in lines
                  if re.match(r"^\s*\d{5,7}", line)]
        if not parsed or any(item is None for item in parsed):
            return False
        columns = list(zip(*(item[1] for item in parsed)))
        ranges = ((0, 60), (-90, 70), (-90, 70), (0, 2000),
                  (-100, 70), (0, 100), (0, 100))
        for values, (lower, upper) in zip(columns, ranges):
            observed = [value for value in values if not math.isclose(value, -99.0)]
            if any(value < lower or value > upper for value in observed):
                return False
        for tmax, tmin in zip(columns[1], columns[2]):
            if (not math.isclose(tmax, -99.0) and not math.isclose(tmin, -99.0)
                    and tmax < tmin):
                return False
        column_names = ("SRAD", "TMAX", "TMIN", "RAIN", "TDEW", "RH2M", "WIND")
        required = {str(name).upper() for name in (required_columns or ())}
        for name, values in zip(column_names, columns):
            if name in required and any(math.isclose(value, -99.0) for value in values):
                return False
        dates = [_wth_code_to_date(item[0]) for item in parsed]
        if any(value is None for value in dates) or len(set(dates)) != len(dates):
            return False
        if any(current - previous != timedelta(days=1)
               for previous, current in zip(dates, dates[1:])):
            return False
        return end_year is None or dates[-1].year >= int(end_year)
    except (OSError, UnicodeError, ValueError):
        return False

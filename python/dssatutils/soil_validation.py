"""Fixed-column soil preflight and offline rebuilding of derived USDA files."""
from __future__ import annotations

import math
from pathlib import Path
import re


def _soil_lines_issue(lines):
    if not lines:
        return "SOIL.SOL is empty or unreadable"
    errors = [s.strip() for s in lines if s.startswith(("*SOIL ERROR", "Source missing", "No Soil ID"))]
    if errors:
        return " | ".join(errors)
    headers = [i for i, s in enumerate(lines) if re.match(r"^@\s+SLB\b", s)]
    if not headers:
        return "SOIL.SOL has no @ SLB layer table"
    for h in headers:
        fields = list(re.finditer(r"\S+", lines[h]))[1:]
        ends = [f.end() for f in fields]
        starts = [0] + [end + 1 for end in ends[:-1]]
        depths = []
        for row in lines[h + 1:]:
            if row.startswith(("@", "*")):
                break
            if not row.strip() or row.lstrip().startswith("!"):
                continue
            tokens = row.split()
            try:
                depth = float(row[:ends[0]])
                # A whitespace parser misses the historic one-space row shift:
                # intended SLB=200 becomes 20 in DSSAT's six-column field.
                if depth != float(tokens[0]) or depth <= 0 or not depth.is_integer():
                    raise ValueError
                for field, start, end in zip(fields, starts, ends):
                    if field.group() == "SLMH":
                        continue  # horizon names may be text
                    value = float(row[start:end])
                    if not math.isfinite(value):
                        raise ValueError
                if any(len(row) > end and not row[end].isspace() for end in ends[:-1]):
                    raise ValueError
            except (ValueError, IndexError):
                return "SOIL.SOL has invalid fixed-width layer fields; regenerate the derived SOL with the corrected writer"
            depths.append(int(depth))
        if not depths:
            return "SOIL.SOL has no parseable SLB layer depths"
        if len(depths) > 19:
            return f"SOIL.SOL has {len(depths)} layers; DSSAT accepts at most 19"
        if any(b <= a for a, b in zip(depths, depths[1:])):
            return "SOIL.SOL layer depths are not strictly increasing: " + ",".join(map(str, depths))
    return None


def soil_file_issue(path):
    """Return None for valid layer columns, otherwise a diagnostic (no mutation).

    Checks formatting, not agronomic plausibility. DSSAT's -99 sentinel is valid.
    """
    path = Path(path)
    if not path.is_file():
        return "SOIL.SOL is missing"
    try:
        return _soil_lines_issue(path.read_text(encoding="utf-8").splitlines())
    except (OSError, UnicodeError):
        return "SOIL.SOL is empty or unreadable"


def rebuild_soil_files_from_mapping(mapping_csv, output_dir, soil_source):
    """Reformat SSURGO/GNATSGO mappings into a NEW directory, without downloads.

    Preserves stored hydraulic/property values; does not rerun pedotransfer
    functions. Refuses an existing destination, including an empty directory.
    Returns a DataFrame with ID and path columns.
    """
    import pandas as pd
    from . import soil_ssurgo, soil_gnatsgo

    source = soil_source.upper()
    if source not in ("SSURGO", "GNATSGO"):
        raise ValueError("soil_source must be SSURGO or GNATSGO")
    destination = Path(output_dir)
    if destination.exists():
        raise ValueError("output_dir must be a new directory; existing caches are never overwritten")
    data = pd.read_csv(mapping_csv, dtype={"ID": str})
    required = {"ID", "latitude", "longitude", "SLLL", "SDUL", "SSAT", "bulk_density", "om_pct", "clay_pct", "silt_pct"}
    if data.empty or not required.issubset(data.columns):
        raise ValueError("soil mapping is empty or lacks required profile columns")
    if data["ID"].isna().any() or not data["ID"].str.fullmatch(r"[A-Za-z0-9_-]{1,10}").all():
        raise ValueError("soil mapping contains invalid IDs")
    if not {"depth_range", "depth_bottom"}.intersection(data.columns):
        raise ValueError("soil mapping lacks layer depths")
    if "depth_bottom" not in data:
        data["depth_bottom"] = pd.to_numeric(data["depth_range"].str.extract(r"-([0-9.]+)cm$", expand=False), errors="raise")
    destination.mkdir(parents=True)
    writer = soil_ssurgo._write_sol if source == "SSURGO" else soil_gnatsgo._write_sol
    records = []
    for soil_id, profile in data.groupby("ID", sort=False):
        writer(profile, str(destination))
        path = destination / f"{soil_id}.SOL"
        issue = soil_file_issue(path)
        if issue:
            raise ValueError(f"{soil_id}: {issue}")
        records.append({"ID": soil_id, "path": str(path)})
    return pd.DataFrame(records, columns=["ID", "path"])

"""Shared configuration helpers for dssatutils."""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml


def _config_candidates() -> list[Path]:
    here = Path(__file__).resolve()
    raw = [
        str(here.parent / "config.yml"),
        str(here.parents[2] / "config.yml"),
        str(Path.cwd() / "config.yml"),
        str(Path.cwd() / "config.yaml"),
        os.environ.get("DSSATUTILS_CONFIG", ""),
    ]
    seen: set[str] = set()
    out: list[Path] = []
    for item in raw:
        if not item:
            continue
        p = Path(item).expanduser()
        key = str(p)
        if key in seen:
            continue
        seen.add(key)
        out.append(p)
    return out


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(merged.get(key), dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


@lru_cache(maxsize=1)
def load_config() -> dict[str, Any]:
    """Load package defaults, then overlay working-directory/user config."""
    config: dict[str, Any] = {}
    for candidate in _config_candidates():
        if not candidate.exists():
            continue
        with candidate.open("r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
        if isinstance(data, dict):
            config = _deep_merge(config, data)
    return config


def get_config_value(path: str, default: Any = None) -> Any:
    """Return a dotted-path config value, or *default* when missing."""
    value: Any = load_config()
    for part in path.split("."):
        if not isinstance(value, dict) or part not in value:
            return default
        value = value[part]
    return default if value is None else value


def get_config_bool(path: str, default: bool = False) -> bool:
    value = get_config_value(path, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", "y"}
    return bool(default)


def get_config_number(path: str, default: float) -> float:
    value = get_config_value(path, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)

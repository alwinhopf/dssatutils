"""Credential helpers for optional remote providers used by dssatutils."""

from __future__ import annotations

import getpass
import os
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional

from .config import get_config_value

CDS_DEFAULT_URL = get_config_value(
    "cds.default_url",
    "https://cds.climate.copernicus.eu/api",
)


def cdsapirc_candidates() -> Iterable[Path]:
    """Yield candidate locations for a cdsapi-compatible ``.cdsapirc`` file."""
    raw_paths = [
        os.environ.get("CDSAPI_RC", ""),
        str(Path.home() / ".cdsapirc"),
        os.path.join(os.environ.get("USERPROFILE", ""), ".cdsapirc"),
        os.path.join(os.environ.get("HOME", ""), ".cdsapirc"),
    ]
    seen = set()
    for raw in raw_paths:
        if not raw:
            continue
        path = Path(raw).expanduser()
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        yield path


def read_cdsapirc(path: Optional[os.PathLike | str] = None) -> Optional[Dict[str, str]]:
    """Read CDS API credentials from ``.cdsapirc`` if present."""
    paths = [Path(path).expanduser()] if path else list(cdsapirc_candidates())
    for candidate in paths:
        if not candidate.exists():
            continue
        values: Dict[str, str] = {}
        with candidate.open("r", encoding="utf-8") as fh:
            for raw in fh:
                if ":" not in raw:
                    continue
                key, value = raw.split(":", 1)
                values[key.strip().lower()] = value.strip()
        token = values.get("key", "")
        if token:
            return {
                "url": values.get("url", CDS_DEFAULT_URL) or CDS_DEFAULT_URL,
                "key": token,
                "path": str(candidate),
            }
    return None


def setup_cds_credentials(
    token: Optional[str] = None,
    url: str = CDS_DEFAULT_URL,
    rc_path: Optional[os.PathLike | str] = None,
    overwrite: bool = False,
    prompt: bool = True,
) -> Dict[str, str]:
    """Configure Copernicus CDS credentials for Python CDS-backed sources.

    The helper accepts an explicit Personal Access Token, uses ``CDSAPI_KEY`` /
    ``CDSAPI_URL``, imports an existing ``.cdsapirc``, or prompts in an
    interactive terminal. It writes a cdsapi-compatible ``.cdsapirc`` when a new
    token is supplied and returns metadata without printing the token.
    """
    token = token or os.environ.get("CDSAPI_KEY", "")
    env_url = os.environ.get("CDSAPI_URL", "")
    if env_url and url == CDS_DEFAULT_URL:
        url = env_url

    existing = None
    if not token and not overwrite:
        existing = read_cdsapirc(rc_path)
        if existing:
            token = existing["key"]
            url = existing.get("url", url) or url
            rc_path = existing["path"]

    if not token and prompt and sys.stdin is not None and sys.stdin.isatty():
        print(
            "Copernicus CDS credentials are required. Create a Personal Access "
            "Token at https://cds.climate.copernicus.eu/how-to-api"
        )
        token = getpass.getpass("Enter Copernicus CDS Personal Access Token: ").strip()

    if not token:
        checked = ", ".join(str(p) for p in cdsapirc_candidates())
        raise RuntimeError(
            "Copernicus CDS credentials were not found. Set CDSAPI_KEY/CDSAPI_URL, "
            "create ~/.cdsapirc, or call setup_cds_credentials(token='<PAT>'). "
            f"Checked: {checked}"
        )

    if rc_path is None:
        rc_path = os.environ.get("CDSAPI_RC") or str(Path.home() / ".cdsapirc")
    rc = Path(rc_path).expanduser()

    if overwrite or not rc.exists():
        rc.parent.mkdir(parents=True, exist_ok=True)
        rc.write_text(f"url: {url}\nkey: {token}\n", encoding="utf-8")
        try:
            rc.chmod(0o600)
        except OSError:
            pass

    os.environ["CDSAPI_URL"] = url
    os.environ["CDSAPI_KEY"] = token
    os.environ["CDSAPI_RC"] = str(rc)
    return {"url": url, "path": str(rc), "has_key": "true"}


def era5land_set_cds_key(
    token: str,
    rc_path: Optional[os.PathLike | str] = None,
    overwrite: bool = True,
) -> Dict[str, str]:
    """Backwards-compatible alias for older ERA5-Land setup scripts."""
    return setup_cds_credentials(
        token=token,
        rc_path=rc_path,
        overwrite=overwrite,
        prompt=False,
    )


def ensure_cds_credentials(prompt: bool = True) -> Dict[str, str]:
    """Ensure CDS credentials are available, prompting only in interactive runs."""
    env_key = os.environ.get("CDSAPI_KEY", "")
    if env_key:
        url = os.environ.get("CDSAPI_URL", CDS_DEFAULT_URL) or CDS_DEFAULT_URL
        return {"url": url, "path": os.environ.get("CDSAPI_RC", ""), "has_key": "true"}

    rc = read_cdsapirc()
    if rc:
        os.environ.setdefault("CDSAPI_URL", rc.get("url", CDS_DEFAULT_URL))
        os.environ.setdefault("CDSAPI_KEY", rc["key"])
        os.environ.setdefault("CDSAPI_RC", rc["path"])
        return {"url": rc.get("url", CDS_DEFAULT_URL), "path": rc["path"], "has_key": "true"}

    return setup_cds_credentials(prompt=prompt)


def make_cds_client(cdsapi, prompt: bool = True):
    """Create a ``cdsapi.Client`` using env vars, ``.cdsapirc``, or a prompt."""
    cfg = ensure_cds_credentials(prompt=prompt)
    kwargs = {"key": os.environ["CDSAPI_KEY"]}
    url = cfg.get("url") or os.environ.get("CDSAPI_URL")
    if url:
        kwargs["url"] = url
    return cdsapi.Client(**kwargs)

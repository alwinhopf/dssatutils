# File: soil_ssurgo_alderman.py
#
# SSURGO querying and processing implementation using full-profile logic
# (dominant component, measured tension fallback, Saxton & Rawls PTFs).
# Mirrors the R soil_ssurgo_alderman.R.

import math
import os
import time
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import requests

try:
    from tqdm import tqdm
    _HAS_TQDM = True
except ImportError:
    _HAS_TQDM = False

_SDA_URL = "https://SDMDataAccess.sc.egov.usda.gov/Tabular/post.rest"

# ---- SDA REST helpers ----------------------------------------------------

def _sda_query(sql: str, max_retries: int = 3, delay: float = 5.0) -> Optional[pd.DataFrame]:
    """POST a SQL query to SDA and return a DataFrame, or None on failure."""
    for attempt in range(max_retries):
        try:
            r = requests.post(
                _SDA_URL,
                data={"query": sql, "format": "json+columnname"},
                timeout=120,
            )
            r.raise_for_status()
            payload = r.json()
            table = payload.get("Table")
            if not table:
                return None
            rows = table[1:]  # first row is headers
            cols = table[0]
            return pd.DataFrame(rows, columns=cols)
        except Exception as exc:
            if attempt < max_retries - 1:
                time.sleep(delay)
            else:
                warnings.warn(f"SDA query failed: {exc}")
                return None


def _sda_spatial_mukeys(lat: float, lon: float, max_retries: int = 3, delay: float = 5.0, STATSGO: bool = False) -> Optional[list]:
    """Return list of mukeys intersecting the given WGS84 point."""
    wkt = f"POINT({lon} {lat})"
    source_condition = "AND lgd.areaname = 'United States'" if STATSGO else "AND lgd.areaname != 'United States'"
    sql = f"""
    SELECT DISTINCT mu.mukey FROM mapunit mu
    INNER JOIN legend lgd ON mu.lkey = lgd.lkey
    WHERE mu.mukey IN (SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('{wkt}'))
    {source_condition}
    """
    df = _sda_query(sql, max_retries=max_retries, delay=delay)
    if df is None or df.empty:
        return None
    return df["mukey"].dropna().tolist()


def _format_in(values) -> str:
    """Format a list as SQL IN clause: ('a','b','c')."""
    cleaned = list({str(v) for v in values if v is not None})
    if not cleaned:
        return "('')"
    inner = ",".join(f"'{v}'" for v in cleaned)
    return f"({inner})"


def _append_log_line(log_file, level="INFO", context="SSURGO", msg="", point_id=None):
    if not log_file:
        return
    import datetime
    pid = os.getpid()
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    id_part = f" [ID={point_id}]" if point_id else ""
    line = f"[{timestamp}] [{level}] [{context}]{id_part} [PID={pid}] {msg}\n"
    try:
        with open(log_file, "a") as fh:
            fh.write(line)
    except:
        pass


def _soil_helper_log(log_file=None, level="INFO", context="SSURGO", msg="", point_id=None):
    _append_log_line(log_file, level, context, msg, point_id)


# ---- Pedotransfer functions ----------------------------------------------

def _ptf_saxton_slll(silt, clay, soc, bd, cf) -> float:
    sand = 1.0 - silt / 100.0 - clay / 100.0
    c_dec = clay / 100.0
    som = soc * 1.72
    theta_1500t = (-0.024 * sand + 0.487 * c_dec + 0.006 * som + 0.005 * sand * som
                   - 0.013 * c_dec * som + 0.068 * sand * c_dec + 0.031)
    theta_1500 = theta_1500t + 0.14 * theta_1500t - 0.02
    a = bd / 2.65
    denom = 1.0 - cf * (1.0 - a)
    if denom == 0:
        denom = 1e-4
    Rv = (a * cf) / denom
    return float(theta_1500 * (1.0 - Rv))


def _ptf_saxton_sdul(silt, clay, soc, bd, cf) -> float:
    sand = 1.0 - silt / 100.0 - clay / 100.0
    c_dec = clay / 100.0
    som = soc * 1.72
    theta_33t = (-0.251 * sand + 0.195 * c_dec + 0.011 * som + 0.006 * sand * som
                 - 0.013 * c_dec * som + 0.452 * sand * c_dec + 0.299)
    theta_33 = theta_33t + 1.283 * (theta_33t ** 2) - 0.374 * theta_33t - 0.015
    a = bd / 2.65
    denom = 1.0 - cf * (1.0 - a)
    if denom == 0:
        denom = 1e-4
    Rv = (a * cf) / denom
    return float(theta_33 * (1.0 - Rv))


def _ptf_saxton_ssat(silt, clay, soc, bd, cf) -> float:
    sand = 1.0 - silt / 100.0 - clay / 100.0
    c_dec = clay / 100.0
    som = soc * 1.72
    theta_33t = (-0.251 * sand + 0.195 * c_dec + 0.011 * som + 0.006 * sand * som
                 - 0.013 * c_dec * som + 0.452 * sand * c_dec + 0.299)
    theta_33 = theta_33t + 1.283 * (theta_33t ** 2) - 0.374 * theta_33t - 0.015
    theta_S33t = (0.278 * sand + 0.034 * c_dec + 0.022 * som - 0.018 * sand * som
                  - 0.027 * c_dec * som - 0.584 * sand * c_dec + 0.078)
    theta_S33 = theta_S33t + 0.636 * theta_S33t - 0.107
    theta_S = theta_33 + theta_S33 - 0.097 * sand + 0.043
    a = bd / 2.65
    denom = 1.0 - cf * (1.0 - a)
    if denom == 0:
        denom = 1e-4
    Rv = (a * cf) / denom
    return float(theta_S * (1.0 - Rv))


def _ptf_saxton_ssks(theta_s, theta_33, theta_1500, cf, bd) -> float:
    theta_s = max(theta_s, 0.02)
    theta_33 = max(theta_33, 0.01)
    theta_1500 = max(theta_1500, 0.005)
    if theta_33 <= 0 or theta_1500 <= 0 or theta_33 <= theta_1500:
        lambda_val = 0.5
    else:
        lambda_val = (math.log(theta_33) - math.log(theta_1500)) / (math.log(1500.0) - math.log(33.0))

    diff_sat_33 = max(theta_s - theta_33, 1e-4)
    Ks = 1930.0 * (diff_sat_33 ** (3.0 - lambda_val))

    denom = 1.0 - cf * (1.0 - 3.0 * (bd / 2.65) / 2.0)
    if denom == 0:
        denom = 1e-4
    Kb = Ks * (1.0 - cf) / denom / 10.0
    return float(max(Kb, 0.001))


def _ptf_slu1(sat, pwp, depth) -> float:
    sat_crop = []
    pwp_crop = []
    depth_crop = []
    for s, p, d in zip(sat, pwp, depth):
        if d > 15.0:
            sat_crop.append(s)
            pwp_crop.append(p)
            depth_crop.append(15.0)
            break
        else:
            sat_crop.append(s)
            pwp_crop.append(p)
            depth_crop.append(d)

    if not depth_crop:
        return 6.0

    val = 0.0
    prev_d = 0.0
    for s, p, d in zip(sat_crop, pwp_crop, depth_crop):
        thick = d - prev_d
        if thick > 0 and not math.isnan(s) and not math.isnan(p):
            val += ((s - p) / 2.0) * thick * 10.0
        prev_d = d
    return float(val)


def _ptf_nrcs_hsg(ksat, depth) -> int:
    ksat_50 = 999.0
    ksat_100 = 999.0
    sl_depth = 0.0
    for k, d in zip(ksat, depth):
        if math.isnan(k):
            continue
        sl_depth = max(sl_depth, d)

    for idx, d in enumerate(depth):
        if d > 50.0:
            ksat_50 = min(ksat[:idx+1])
            break
    else:
        if ksat:
            ksat_50 = min(ksat)

    for idx, d in enumerate(depth):
        if d > 100.0:
            ksat_100 = min(ksat[:idx+1])
            break
    else:
        if ksat:
            ksat_100 = min(ksat)

    hsg = 2
    if ksat_50 == 999.0: ksat_50 = float('nan')
    if ksat_100 == 999.0: ksat_100 = float('nan')

    if not math.isnan(ksat_50) and not math.isnan(ksat_100) and sl_depth > 0:
        if (ksat_50 > 40.0 * 0.36 and 50.0 <= sl_depth <= 100.0) or (ksat_100 > 10.0 * 0.36 and sl_depth > 100.0):
            hsg = 1
        elif (40.0 * 0.36 >= ksat_50 > 10.0 * 0.36 and 50.0 <= sl_depth <= 100.0) or (10.0 * 0.36 >= ksat_100 > 4.0 * 0.36 and sl_depth > 100.0):
            hsg = 2
        elif (10.0 * 0.36 >= ksat_50 > 1.0 * 0.36 and 50.0 <= sl_depth <= 100.0) or (4.0 * 0.36 >= ksat_100 > 0.4 * 0.36 and sl_depth > 100.0):
            hsg = 3
        elif (ksat_50 <= 1.0 * 0.36 and 50.0 <= sl_depth <= 100.0) or sl_depth < 50.0 or (ksat_100 <= 0.4 * 0.36 and sl_depth > 100.0):
            hsg = 4
    return hsg


def _ptf_curve_number(slope, hsg, ksat=None, depth=None) -> float:
    if math.isnan(hsg):
        hsg = _ptf_nrcs_hsg(ksat, depth)

    cn = 73
    if hsg == 1:
        if 0 <= slope <= 2: cn = 61
        elif 2 < slope <= 5: cn = 64
        elif 5 < slope <= 10: cn = 68
        elif slope > 10: cn = 71
    elif hsg == 2:
        if 0 <= slope <= 2: cn = 73
        elif 2 < slope <= 5: cn = 76
        elif 5 < slope <= 10: cn = 80
        elif slope > 10: cn = 83
    elif hsg == 3:
        if 0 <= slope <= 2: cn = 81
        elif 2 < slope <= 5: cn = 84
        elif 5 < slope <= 10: cn = 88
        elif slope > 10: cn = 91
    elif hsg == 4:
        if 0 <= slope <= 2: cn = 84
        elif 2 < slope <= 5: cn = 87
        elif 5 < slope <= 10: cn = 91
        elif slope > 10: cn = 94
    return float(cn)


def _sldr_from_drainage(drainage) -> float:
    d = str(drainage).strip()
    mapping = {
        "Excessively drained": 0.85,
        "Somewhat excessively drained": 0.75,
        "Well drained": 0.60,
        "Moderately well drained": 0.40,
        "Somewhat poorly drained": 0.25,
        "Poorly drained": 0.05,
        "Very poorly drained": 0.01
    }
    return mapping.get(d, 0.60)


# ---- Helper functions for horizon conversion -----------------------------

def _coalesce_num(*args) -> Optional[float]:
    for a in args:
        if a is not None and not (isinstance(a, float) and math.isnan(a)) and a != "":
            try:
                return float(a)
            except ValueError:
                pass
    return None


def _clip01(x) -> Optional[float]:
    if x is None or math.isnan(x):
        return None
    return min(max(x, 0.001), 0.95)


def _calc_sbdm(row) -> Optional[float]:
    return _coalesce_num(row.get("dbtenthbar_r"), row.get("dbthirdbar_r"), row.get("dbovendry_r"))


def _is_coarse_texture(sand, silt, clay) -> bool:
    if math.isnan(sand) or math.isnan(silt) or math.isnan(clay):
        return False
    return (
        (sand >= 85 and (silt + 1.5 * clay) <= 15) or
        (sand >= 85 and sand < 90 and (silt + 1.5 * clay) >= 15) or
        (sand >= 70 and sand < 85 and (silt + 2 * clay) <= 30) or
        (clay <= 20 and sand >= 52 and (silt + 2 * clay) > 30) or
        (clay < 7 and silt < 50 and sand > 43 and sand < 52)
    )


def _calc_sdul_measured(row) -> Optional[float]:
    clay = _coalesce_num(row.get("claytotal_r"))
    sand = _coalesce_num(row.get("sandtotal_r"))
    silt = _coalesce_num(row.get("silttotal_r"))
    coarse = _is_coarse_texture(sand or 0, silt or 0, clay or 0)
    wthird = _coalesce_num(row.get("wthirdbar_r"))
    if wthird is not None:
        wthird /= 100.0
    else:
        return None
    wtenth = _coalesce_num(row.get("wtenthbar_r"))
    if wtenth is not None:
        wtenth /= 100.0
    wtenth = _coalesce_num(wtenth, wthird)
    return wtenth if coarse else wthird


def _calc_ssat_measured(row) -> Optional[float]:
    partdensity = _coalesce_num(row.get("partdensity"), 2.65)
    wtenthbar = _coalesce_num(row.get("wtenthbar_r"))
    wthirdbar = _coalesce_num(row.get("wthirdbar_r"))
    dbtenthbar = _coalesce_num(row.get("dbtenthbar_r")) if wtenthbar is not None else None
    dbthirdbar = _coalesce_num(row.get("dbthirdbar_r")) if wthirdbar is not None else None

    ssat_tenth = None
    if dbtenthbar is not None:
        ssat_tenth = 0.95 * (1.0 - dbtenthbar / partdensity)

    ssat_third = None
    if dbthirdbar is not None:
        ssat_third = 0.95 * (1.0 - dbthirdbar / partdensity)

    wsatiated = _coalesce_num(row.get("wsatiated_r"))
    ssat_wsat = wsatiated / 100.0 if wsatiated is not None else None

    dbovendry = _coalesce_num(row.get("dbovendry_r"))
    ssat_dry = None
    if dbovendry is not None:
        ssat_dry = 0.95 * (1.0 - dbovendry / partdensity)

    return _coalesce_num(ssat_wsat, ssat_tenth, ssat_third, ssat_dry)


# ---- Dominant component profile builder ----------------------------------

def _build_dssat_profile_from_component(component_row, horizon_tbl: pd.DataFrame, point_id: str, lat: float, lon: float, log_file=None, standardize_layers=False) -> Optional[dict]:
    _soil_helper_log(log_file, "INFO", "SSURGO_DOMINANT",
                     f"Building dominant profile from cokey={component_row.get('cokey')} with {len(horizon_tbl)} raw horizon row(s)",
                     point_id)

    if standardize_layers:
        hz_depb_numeric = pd.to_numeric(horizon_tbl["hzdepb_r"], errors="coerce")
        bedrock_depth = hz_depb_numeric.max() if not hz_depb_numeric.dropna().empty else 200.0
        if not math.isfinite(bedrock_depth) or pd.isna(bedrock_depth) or bedrock_depth <= 0:
            bedrock_depth = 200.0

        std_horizon_tbl = _aggregate_horizons_to_standard_layers(horizon_tbl, bedrock_depth)
        if std_horizon_tbl is None or std_horizon_tbl.empty:
            _soil_helper_log(log_file, "ERROR", "SSURGO_DOMINANT", "Layer standardization aggregated 0 layers", point_id)
            return None
        horizon_tbl = std_horizon_tbl

    hz_df = horizon_tbl.copy()
    num_cols = ["hzdept_r", "hzdepb_r", "dbovendry_r", "dbtenthbar_r", "dbthirdbar_r",
                "dbfifteenbar_r", "wsatiated_r", "wtenthbar_r", "wthirdbar_r",
                "partdensity", "ksat_r", "wfifteenbar_r", "sandtotal_r", "claytotal_r",
                "silttotal_r", "om_r", "fragvol_r"]
    for col in num_cols:
        if col in hz_df.columns:
            hz_df[col] = pd.to_numeric(hz_df[col], errors="coerce")

    hz_df = hz_df.sort_values(by="hzdepb_r")

    def agg_frag(x):
        non_na = x.dropna()
        if non_na.empty:
            return np.nan
        return non_na.sum()

    agg_dict = {
        "hzdept_r": "min",
        "fragvol_r": agg_frag
    }
    for col in num_cols:
        if col not in ["hzdept_r", "fragvol_r", "hzdepb_r"]:
            agg_dict[col] = "mean"

    # Group by hzdepb_r, hzname, cokey
    # fill missing hzname or cokey to avoid dropping them
    hz_df["hzname"] = hz_df["hzname"].fillna("-99")
    hz_df["cokey"] = hz_df["cokey"].fillna("-99")

    grouped = hz_df.groupby(["hzdepb_r", "hzname", "cokey"], as_index=False).agg(agg_dict)
    grouped = grouped.sort_values(by="hzdepb_r").reset_index(drop=True)
    grouped = grouped.replace([np.inf, -np.inf], np.nan)

    layers_list = []
    for idx, row in grouped.iterrows():
        hzname = str(row.get("hzname") or "")
        bedrock = any(c in hzname for c in ["r", "R"])

        fragvol_raw = row.get("fragvol_r")
        if fragvol_raw is not None and not math.isnan(fragvol_raw):
            fragvol_raw = min(fragvol_raw, 99.0)

        fragvol_r = _coalesce_num(fragvol_raw, 0.0)
        coarse_fraction = fragvol_r / 100.0
        partdensity = _coalesce_num(row.get("partdensity"), 2.65)

        SBDM = _calc_sbdm(row)
        SSAT_raw = _calc_ssat_measured(row)
        SDUL_raw = _calc_sdul_measured(row)

        wfifteen = _coalesce_num(row.get("wfifteenbar_r"))
        SLLL_raw = wfifteen / 100.0 if wfifteen is not None else None

        om = _coalesce_num(row.get("om_r"), 0.0)
        soc = om / 1.724

        clay = _coalesce_num(row.get("claytotal_r"), 0.0)
        silt = _coalesce_num(row.get("silttotal_r"), 0.0)
        sand = 100.0 - clay - silt

        if bedrock:
            SLLL_ptf = np.nan
            SDUL_ptf = np.nan
            SSAT_ptf = np.nan
        else:
            bd_for_ptf = _coalesce_num(SBDM, row.get("dbthirdbar_r"), row.get("dbovendry_r"), 1.4)
            SLLL_ptf = _ptf_saxton_slll(silt, clay, soc, bd_for_ptf, coarse_fraction)
            SDUL_ptf = _ptf_saxton_sdul(silt, clay, soc, bd_for_ptf, coarse_fraction)
            SSAT_ptf = _ptf_saxton_ssat(silt, clay, soc, bd_for_ptf, coarse_fraction)

        # Floor wilting point at 0.02 (consistent with soil_ssurgo / soil_gnatsgo):
        # measured wfifteenbar_r or the Saxton-Rawls fallback can yield SLLL ~0 on
        # sandy soils, which SIGFPEs DSSAT. The sdul = max(sdul, slll + 0.04) step
        # in the fill loop keeps a usable plant-available-water gap above the
        # floored LL (a gap < ~0.04 still SIGFPEs DSSAT's water balance mid-season).
        SLLL = max(_clip01(_coalesce_num(SLLL_raw, SLLL_ptf)), 0.02)
        SDUL = _clip01(_coalesce_num(SDUL_raw, SDUL_ptf))
        SSAT = _clip01(_coalesce_num(SSAT_raw, SSAT_ptf))

        layers_list.append({
            "hzdepb_r": row.get("hzdepb_r"),
            "hzdept_r": row.get("hzdept_r"),
            "hzname": hzname,
            "bedrock": bedrock,
            "fragvol_raw": fragvol_raw,
            "coarse_fraction": coarse_fraction,
            "partdensity": partdensity,
            "soc": soc,
            "SBDM": SBDM,
            "SLLL": SLLL,
            "SDUL": SDUL,
            "SSAT": SSAT,
            "clay": clay,
            "silt": silt,
            "sand": sand,
            "ksat_r": row.get("ksat_r"),
            "dbthirdbar_r": row.get("dbthirdbar_r"),
            "dbovendry_r": row.get("dbovendry_r")
        })

    layers_df = pd.DataFrame(layers_list)
    if layers_df.empty:
        return None

    cols_to_fill = ["SDUL", "SLLL", "SSAT", "silt", "clay", "SBDM", "soc"]
    layers_df[cols_to_fill] = layers_df[cols_to_fill].ffill()

    processed_layers = []
    for idx, row in layers_df.iterrows():
        slll = row["SLLL"]
        sdul = row["SDUL"]
        ssat = row["SSAT"]

        if slll is not None and not math.isnan(slll):
            if sdul is not None and not math.isnan(sdul):
                sdul = max(sdul, slll + 0.04)
            if ssat is not None and not math.isnan(ssat):
                ssat = max(ssat, (sdul or slll) + 0.01)

        bedrock = row["bedrock"]
        fragvol_raw = row["fragvol_raw"]

        if bedrock and (fragvol_raw is None or math.isnan(fragvol_raw)):
            slcf = 99.0
        else:
            slcf = _coalesce_num(fragvol_raw, 0.0)

        srgf = max(0.01, 1.0 - slcf / 100.0) if bedrock else 1.0

        ksat_r = _coalesce_num(row["ksat_r"])
        is_ksat_low = ksat_r is not None and ksat_r < 0.001 / 60 / 60 * 10000
        if (bedrock and ksat_r is None) or is_ksat_low:
            ssks = 0.001
        else:
            cf = row["coarse_fraction"]
            bd = _coalesce_num(row["SBDM"], 1.4)
            ssks_ptf = _ptf_saxton_ssks(ssat or 0.5, sdul or 0.3, slll or 0.15, cf, bd)
            ssks = _coalesce_num(ksat_r * 0.36 if ksat_r is not None else None, ssks_ptf)

        slb = int(round(row["hzdepb_r"]))

        processed_layers.append({
            "SLB": slb,
            "SLMH": str(row["hzname"]) if row["hzname"] != "" and row["hzname"] is not None else "-99",
            "SLLL": slll,
            "SDUL": sdul,
            "SSAT": ssat,
            "SRGF": srgf,
            "SSKS": ssks,
            "SBDM": _coalesce_num(row["SBDM"], row["dbthirdbar_r"], row["dbovendry_r"], 1.4),
            "SLOC": _coalesce_num(row["soc"], 0.0),
            "SLCL": _coalesce_num(row["clay"], 0.0),
            "SLSI": _coalesce_num(row["silt"], 0.0),
            "SLCF": slcf,
            "SLNI": float('nan'),
            "SLHW": float('nan'),
            "SLHB": float('nan'),
            "SCEC": float('nan'),
            "SADC": float('nan')
        })

    final_layers_df = pd.DataFrame(processed_layers)
    final_layers_df = final_layers_df[final_layers_df["SLB"] > 0].drop_duplicates(subset="SLB").sort_values(by="SLB").reset_index(drop=True)
    if final_layers_df.empty:
        return None

    if final_layers_df["SLLL"].isna().all() or final_layers_df["SDUL"].isna().all() or final_layers_df["SSAT"].isna().all():
        _soil_helper_log(log_file, "ERROR", "SSURGO_DOMINANT", "Dominant component profile has all NA for hydraulic properties", point_id)
        return None

    _soil_helper_log(log_file, "INFO", "SSURGO_DOMINANT",
                     f"Dominant profile retained {len(final_layers_df)} DSSAT layer(s): "
                     f"{','.join(str(int(b)) for b in final_layers_df['SLB'])}",
                     point_id)

    slu1 = _ptf_slu1(final_layers_df["SSAT"].tolist(), final_layers_df["SLLL"].tolist(), final_layers_df["SLB"].tolist())
    if not math.isfinite(slu1) or math.isnan(slu1):
        slu1 = 6.0

    albedo = _coalesce_num(component_row.get("albedodry_r"), 0.13)
    if albedo <= 0.0:
        albedo = 0.13
    drainage = component_row.get("drainage")
    sldr = _sldr_from_drainage(drainage)

    slope = _coalesce_num(component_row.get("slope_r"), 0.0)
    hydgrp = str(component_row.get("hydgrp") or "")
    hsg_code = float('nan')
    if hydgrp and hydgrp[0].upper() in ["A", "B", "C", "D"]:
        hsg_code = float(["A", "B", "C", "D"].index(hydgrp[0].upper()) + 1)

    slro = _ptf_curve_number(
        slope, hsg_code,
        final_layers_df["SSKS"].tolist(), final_layers_df["SLB"].tolist()
    )

    metadata = pd.DataFrame([{
        "ID": point_id,
        "SOIL_ID": point_id,
        "mukey": str(component_row.get("mukey")),
        "cokey": str(component_row.get("cokey")),
        "compname": str(component_row.get("compname")),
        "comppct_r": _coalesce_num(component_row.get("comppct_r")),
        "latitude": lat,
        "longitude": lon
    }])

    return {
        "profile_id": point_id,
        "latitude": lat,
        "longitude": lon,
        "site": point_id,
        "country": "USA",
        "scs_family": "",
        "scom": "SC",
        "salb": albedo,
        "slu1": round(slu1, 1),
        "sldr": sldr,
        "slro": slro,
        "slnf": 1.0,
        "slpf": 1.0,
        "smhb": "IB001",
        "smpx": "IB001",
        "smke": "IB001",
        "layers": final_layers_df,
        "metadata": metadata
    }


# ---- Weighted simple fallback profile builder ----------------------------

def _calculate_soil_properties_fallback(props_df: pd.DataFrame, top_depth: float, bottom_depth: float) -> Optional[dict]:
    df = props_df.copy()
    for col in ["hzdept_r", "hzdepb_r", "claytotal_r", "sandtotal_r", "om_r", "dbthirdbar_r", "comppct_r"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df["adj_top"] = df["hzdept_r"].clip(lower=top_depth)
    df["adj_bottom"] = df["hzdepb_r"].clip(upper=bottom_depth)
    df["thickness"] = (df["adj_bottom"] - df["adj_top"]).clip(lower=0)
    df = df[df["thickness"] > 0].copy()
    if df.empty:
        return None

    w = df["thickness"] * df["comppct_r"]
    total_w = w.sum()
    if total_w == 0 or np.isnan(total_w):
        return None

    clay_pct = (df["claytotal_r"] * w).sum() / total_w
    sand_pct = (df["sandtotal_r"] * w).sum() / total_w
    om_pct = (df["om_r"] * w).sum() / total_w
    bulk_density = (df["dbthirdbar_r"] * w).sum() / total_w

    return {
        "clay_pct": clay_pct if not np.isnan(clay_pct) else 20.0,
        "sand_pct": sand_pct if not np.isnan(sand_pct) else 40.0,
        "silt_pct": max(0.0, 100.0 - (clay_pct if not np.isnan(clay_pct) else 20.0) - (sand_pct if not np.isnan(sand_pct) else 40.0)),
        "om_pct": om_pct if not np.isnan(om_pct) else 1.0,
        "bulk_density": bulk_density if not np.isnan(bulk_density) else 1.4,
        "depth_range": f"{top_depth}-{bottom_depth}cm"
    }


def _build_simple_fallback_profile(lat: float, lon: float, point_id: str, mukeys: list, comp_df: Optional[pd.DataFrame], log_file=None, STATSGO: bool = False) -> Optional[dict]:
    _soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK", "Starting weighted-layer fallback profile build", point_id)

    if mukeys:
        source_condition = "lgd.areaname = 'United States'" if STATSGO else "lgd.areaname != 'United States'"
        q_filter = f"""
        SELECT mu.mukey FROM mapunit mu
        INNER JOIN legend lgd ON mu.lkey = lgd.lkey
        WHERE mu.mukey IN {_format_in(mukeys)}
        AND {source_condition}
        """
        res_filter = _sda_query(q_filter)
        if res_filter is not None and not res_filter.empty:
            mukeys = res_filter["mukey"].dropna().astype(str).unique().tolist()
        else:
            mukeys = []

    # query bedrock
    q_bedrock = f"SELECT mukey, brockdepmin FROM muaggatt WHERE mukey IN {_format_in(mukeys)}"
    bedrock_data = _sda_query(q_bedrock)
    bedrock_depth = 200.0
    if bedrock_data is not None and not bedrock_data.empty:
        bd_vals = pd.to_numeric(bedrock_data["brockdepmin"], errors="coerce").dropna().tolist()
        if bd_vals:
            bedrock_depth = min(bd_vals)

    if not math.isfinite(bedrock_depth) or bedrock_depth <= 0:
        bedrock_depth = 200.0

    all_layers = [
        ('0-5cm', (0, 5)), ('5-20cm', (5, 20)), ('20-35cm', (20, 35)),
        ('35-50cm', (35, 50)), ('50-65cm', (50, 65)), ('65-80cm', (65, 80)),
        ('80-95cm', (80, 95)), ('95-110cm', (95, 110)), ('110-125cm', (110, 125)),
        ('125-140cm', (125, 140)), ('140-155cm', (140, 155)), ('155-170cm', (155, 170)),
        ('170-185cm', (170, 185)), ('185-200cm', (185, 200))
    ]
    valid_layers = [(name, limits) for name, limits in all_layers if limits[0] < bedrock_depth]
    if valid_layers:
        name, limits = valid_layers[-1]
        if limits[1] > bedrock_depth:
            valid_layers[-1] = (f"{limits[0]}-{int(bedrock_depth)}cm", (limits[0], bedrock_depth))
    else:
        valid_layers = [(f"0-{int(bedrock_depth)}cm", (0, bedrock_depth))]

    _soil_helper_log(log_file, "INFO", "SSURGO_FALLBACK",
                     f"Fallback using bedrock depth {bedrock_depth:.1f} cm and {len(valid_layers)} target layer(s)",
                     point_id)

    q_soil = (
        "SELECT component.mukey, component.cokey, component.comppct_r, "
        "chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.claytotal_r, "
        "chorizon.sandtotal_r, chorizon.om_r, chorizon.dbthirdbar_r "
        "FROM component INNER JOIN chorizon ON component.cokey = chorizon.cokey "
        f"WHERE component.mukey IN {_format_in(mukeys)}"
    )
    props = _sda_query(q_soil)
    if props is None or props.empty:
        _soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback property query returned empty/NULL", point_id)
        return None

    layers_list = []
    for name, (top, bot) in valid_layers:
        cp = _calculate_soil_properties_fallback(props, top, bot)
        if cp:
            cp["depth_range"] = name
            layers_list.append(cp)

    if not layers_list:
        _soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback aggregated 0 layers", point_id)
        return None

    results_df = pd.DataFrame(layers_list)

    # Choose dominant component for metadata and site properties
    if comp_df is not None and not comp_df.empty:
        dom_comp = _choose_dominant_component(comp_df)
    else:
        dom_comp = pd.Series({
            "cokey": None, "mukey": mukeys[0], "compname": "SSURGO",
            "comppct_r": np.nan, "hydgrp": "", "slope_r": np.nan,
            "drainage": "Well drained", "albedodry_r": 0.13
        })

    processed_layers = []
    for idx, row in results_df.iterrows():
        depth_str = str(row["depth_range"]).split("-")[1].replace("cm", "")
        depth_num = float(depth_str)
        soc = row["om_pct"] / 1.724

        bd = _coalesce_num(row["bulk_density"], 1.4)
        slll = _clip01(_ptf_saxton_slll(row["silt_pct"], row["clay_pct"], soc, bd, 0.0))
        sdul = _clip01(_ptf_saxton_sdul(row["silt_pct"], row["clay_pct"], soc, bd, 0.0))
        ssat = _clip01(_ptf_saxton_ssat(row["silt_pct"], row["clay_pct"], soc, bd, 0.0))

        if slll is not None and not math.isnan(slll):
            if sdul is not None and not math.isnan(sdul):
                sdul = max(sdul, slll + 0.04)
            if ssat is not None and not math.isnan(ssat):
                ssat = max(ssat, (sdul or slll) + 0.01)

        ssks = _ptf_saxton_ssks(ssat or 0.5, sdul or 0.3, slll or 0.15, 0.0, bd)
        slb = int(round(depth_num))

        processed_layers.append({
            "SLB": slb,
            "SLMH": "-99",
            "SLLL": slll,
            "SDUL": sdul,
            "SSAT": ssat,
            "SRGF": 1.0,
            "SSKS": ssks,
            "SBDM": bd,
            "SLOC": soc,
            "SLCL": row["clay_pct"],
            "SLSI": row["silt_pct"],
            "SLCF": 0.0,
            "SLNI": float('nan'),
            "SLHW": float('nan'),
            "SLHB": float('nan'),
            "SCEC": float('nan'),
            "SADC": float('nan')
        })

    layers_df = pd.DataFrame(processed_layers).sort_values(by="SLB").reset_index(drop=True)
    if layers_df.empty:
        return None

    if layers_df["SLLL"].isna().all() or layers_df["SDUL"].isna().all() or layers_df["SSAT"].isna().all():
        _soil_helper_log(log_file, "ERROR", "SSURGO_FALLBACK", "Fallback profile has all NA for hydraulic properties", point_id)
        return None

    slu1 = _ptf_slu1(layers_df["SSAT"].tolist(), layers_df["SLLL"].tolist(), layers_df["SLB"].tolist())
    if not math.isfinite(slu1) or math.isnan(slu1):
        slu1 = 6.0

    albedo = _coalesce_num(dom_comp.get("albedodry_r"), 0.13)
    if albedo <= 0.0:
        albedo = 0.13
    sldr = _sldr_from_drainage(dom_comp.get("drainage", "Well drained"))

    slope = _coalesce_num(dom_comp.get("slope_r"), 0.0)
    hydgrp = str(dom_comp.get("hydgrp") or "")
    hsg_code = float('nan')
    if hydgrp and hydgrp[0].upper() in ["A", "B", "C", "D"]:
        hsg_code = float(["A", "B", "C", "D"].index(hydgrp[0].upper()) + 1)

    slro = _ptf_curve_number(
        slope, hsg_code,
        layers_df["SSKS"].tolist(), layers_df["SLB"].tolist()
    )

    metadata = pd.DataFrame([{
        "ID": point_id,
        "SOIL_ID": point_id,
        "mukey": ";".join(mukeys),
        "cokey": str(dom_comp.get("cokey")) if dom_comp.get("cokey") is not None else "",
        "compname": str(dom_comp.get("compname", "SSURGO")),
        "comppct_r": _coalesce_num(dom_comp.get("comppct_r")),
        "latitude": lat,
        "longitude": lon
    }])

    return {
        "profile_id": point_id,
        "latitude": lat,
        "longitude": lon,
        "site": point_id,
        "country": "USA",
        "scs_family": "",
        "scom": "SC",
        "salb": albedo,
        "slu1": round(slu1, 1),
        "sldr": sldr,
        "slro": slro,
        "slnf": 1.0,
        "slpf": 1.0,
        "smhb": "IB001",
        "smpx": "IB001",
        "smke": "IB001",
        "layers": layers_df,
        "metadata": metadata
    }


# ---- DSSAT output formatting & file writer -------------------------------

def format_dssat_decimal(x, digits=3, width=5) -> str:
    if x is None or math.isnan(x):
        return f"{'-99':>{width}s}"
    out = f"{x:{width}.{digits}f}"
    if out.startswith("0."):
        out = " " + out[1:]
    return out


def _write_dssat_soil_file(profile, output_dir):
    profile_id = profile["profile_id"]
    filename = os.path.join(output_dir, f"{profile_id}.SOL")

    lines = [
        "*SOILS: USA SSURGO Soil Profiles",
        "! Generated from SSURGO database using full-profile logic",
        "",
        f"*{profile_id:<10s} SSURGO        {profile['latitude']:9.3f} {profile['longitude']:9.3f}",
        "@SITE        COUNTRY          LAT     LONG SCS FAMILY",
        f" {profile['site']:<11s} {profile['country']:<10s} {profile['latitude']:9.3f} {profile['longitude']:9.3f} {profile['scs_family']:<20s}",
        "@ SCOM  SALB  SLU1  SLDR  SLRO  SLNF  SLPF  SMHB  SMPX  SMKE",
        f" {profile['scom']:>5s} {profile['salb']:5.2f} {profile['slu1']:5.1f} {profile['sldr']:5.2f} {profile['slro']:5.0f} {profile['slnf']:5.0f} {profile['slpf']:5.0f} {profile['smhb']:>5s} {profile['smpx']:>5s} {profile['smke']:>5s}",
        "@  SLB  SLMH  SLLL  SDUL  SSAT  SRGF  SSKS  SBDM  SLOC  SLCL  SLSI  SLCF  SLNI  SLHW  SLHB  SCEC  SADC"
    ]

    for _, lyr in profile["layers"].iterrows():
        slb = int(lyr["SLB"])
        slmh = str(lyr["SLMH"])[:5]
        slll = format_dssat_decimal(lyr["SLLL"], 3, 5)
        sdul = format_dssat_decimal(lyr["SDUL"], 3, 5)
        ssat = format_dssat_decimal(lyr["SSAT"], 3, 5)
        srgf = lyr["SRGF"] if not math.isnan(lyr["SRGF"]) else 1.0
        ssks = lyr["SSKS"] if not math.isnan(lyr["SSKS"]) else -99.0
        sbdm = lyr["SBDM"] if not math.isnan(lyr["SBDM"]) else -99.0
        sloc = lyr["SLOC"] if not math.isnan(lyr["SLOC"]) else -99.0
        slcl = lyr["SLCL"] if not math.isnan(lyr["SLCL"]) else -99.0
        slsi = lyr["SLSI"] if not math.isnan(lyr["SLSI"]) else -99.0
        slcf = lyr["SLCF"] if not math.isnan(lyr["SLCF"]) else -99.0

        slni = format_dssat_decimal(lyr["SLNI"], 3, 5) if not math.isnan(lyr["SLNI"]) else "  -99"
        slhw = format_dssat_decimal(lyr["SLHW"], 3, 5) if not math.isnan(lyr["SLHW"]) else "  -99"
        slhb = format_dssat_decimal(lyr["SLHB"], 3, 5) if not math.isnan(lyr["SLHB"]) else "  -99"
        scec = format_dssat_decimal(lyr["SCEC"], 3, 5) if not math.isnan(lyr["SCEC"]) else "  -99"
        sadc = format_dssat_decimal(lyr["SADC"], 3, 5) if not math.isnan(lyr["SADC"]) else "  -99"

        # Build layer line matching R sprint format:
        # "%5d %5s %5s %5s %5s %5.2f %5.2f %5.2f %5.2f %5.1f %5.1f %5.0f %5s %5s %5s %5s %5s"
        line = (
            f"{slb:5d} {slmh:>5s} {slll:>5s} {sdul:>5s} {ssat:>5s} "
            f"{srgf:5.2f} {ssks:5.2f} {sbdm:5.2f} {sloc:5.2f} "
            f"{slcl:5.1f} {slsi:5.1f} {slcf:5.0f} "
            f"{slni:>5s} {slhw:>5s} {slhb:>5s} {scec:>5s} {sadc:>5s}"
        )
        lines.append(line)

    lines.append("")
    with open(filename, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---- Per-point worker ----------------------------------------------------

def _process_point(args: dict) -> Optional[pd.DataFrame]:
    """Query SSURGO for one point and write its .SOL file."""
    point_id = args["ID"]
    lat = args["lat"]
    lon = args["lon"]
    output_dir = args["output_dir"]
    log_file = args.get("log_file")
    soil_name = args.get("soil_name")
    STATSGO = args.get("STATSGO", False)
    standardize_layers = args.get("standardize_layers", False)

    _soil_helper_log(log_file, "INFO", "SSURGO_POINT", "Starting SSURGO point processing", point_id)
    profile = None

    if soil_name:
        profile = pull_profile_by_name_alderman(
            soil_name=soil_name,
            lat=lat,
            lon=lon,
            SSURGO=not STATSGO,
            STATSGO=STATSGO,
            standardize_layers=standardize_layers,
            log_file=log_file,
            profile_id=point_id
        )
    else:
        # 1. Spatial query -> mukeys
        mukeys = _sda_spatial_mukeys(lat, lon, STATSGO=STATSGO)
        comp_df = None
        if mukeys:
            # 2. Get component table
            comp_df = _get_component_table(mukeys, STATSGO=STATSGO)
            if comp_df is not None and not comp_df.empty:
                # Choose dominant component
                comp_row = _choose_dominant_component(comp_df)
                cokey = comp_row.get("cokey")

                # Query horizons for cokey
                horizons = _get_component_horizons(cokey)
                if horizons is not None and not horizons.empty:
                    profile = _build_dssat_profile_from_component(
                        comp_row, horizons, point_id, lat, lon, log_file,
                        standardize_layers=standardize_layers
                    )

        # 3. Fallback if profile construction failed
        if profile is None and mukeys:
            _soil_helper_log(log_file, "WARN", "SSURGO_POINT", "Falling back to weighted-layer SSURGO profile logic", point_id)
            profile = _build_simple_fallback_profile(lat, lon, point_id, mukeys, comp_df, log_file, STATSGO=STATSGO)

    if profile is None:
        _soil_helper_log(log_file, "ERROR", "SSURGO_POINT", "No soil profile could be generated for this point", point_id)
        return None

    # Write SOL file
    _write_dssat_soil_file(profile, output_dir)
    _soil_helper_log(log_file, "INFO", "SSURGO_POINT", f"Wrote soil profile to {os.path.join(output_dir, f'{point_id}.SOL')}", point_id)
    return profile["metadata"]


def _get_component_table(mukeys: list, STATSGO: bool = False) -> Optional[pd.DataFrame]:
    source_condition = "legend.areaname = 'United States'" if STATSGO else "legend.areaname != 'United States'"
    sql = f"""
    SELECT component.compname, component.cokey, component.mukey, COALESCE(component.comppct_r,'') AS comppct_r,
           COALESCE(component.hydgrp,'') AS hydgrp, COALESCE(component.slope_r,'') AS slope_r,
           COALESCE(component.drainagecl,'') AS drainage, COALESCE(component.albedodry_r,'') AS albedodry_r
    FROM component
    INNER JOIN mapunit ON component.mukey = mapunit.mukey
    INNER JOIN legend ON mapunit.lkey = legend.lkey
    WHERE component.mukey IN {_format_in(mukeys)}
    AND {source_condition}
    """
    return _sda_query(sql)


def _choose_dominant_component(comp_df: pd.DataFrame) -> pd.Series:
    df = comp_df.copy()
    df["comppct_r"] = pd.to_numeric(df["comppct_r"], errors="coerce")
    df = df.sort_values(by=["comppct_r", "cokey"], ascending=[False, True])
    return df.iloc[0]


def _get_component_horizons(cokey: str) -> Optional[pd.DataFrame]:
    sql = f"""
    SELECT chorizon.hzdept_r, chorizon.hzdepb_r, chorizon.dbovendry_r,
           chorizon.dbtenthbar_r, chorizon.dbthirdbar_r, chorizon.dbfifteenbar_r,
           chorizon.wsatiated_r, chorizon.wtenthbar_r, chorizon.wthirdbar_r,
           chorizon.partdensity, chorizon.ksat_r, chorizon.wfifteenbar_r,
           chorizon.sandtotal_r, chorizon.claytotal_r, chorizon.silttotal_r,
           chorizon.om_r, chorizon.hzname, chfrags.fragvol_r AS fragvol_r,
           chorizon.cokey FROM chorizon LEFT JOIN chfrags ON chfrags.chkey = chorizon.chkey
    WHERE chorizon.cokey IN ('{cokey}') ORDER BY chorizon.hzdepb_r
    """
    return _sda_query(sql)


# ---- Public entry point --------------------------------------------------

def _filter_mukey_by_coord_alderman(mukey_df: pd.DataFrame, lat: float, lon: float, soil_name: str) -> str:
    import geopandas as gpd
    from shapely.geometry import Point
    from shapely.wkt import loads

    # Radii in meters
    radii = [1e3, 1e4, 5e4] + [1e5 * i for i in range(1, 11)] + [1e6 + 1e5 * i for i in range(1, 11)]

    pt = Point(lon, lat)
    pt_gdf = gpd.GeoDataFrame(geometry=[pt], crs="EPSG:4326")
    # Transform to Albers Equal Area (same projection as in R)
    pt_aea = pt_gdf.to_crs("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")

    mukey_subset = None
    mukeys_in = mukey_df["mukey"].dropna().astype(str).unique().tolist()

    for radius in radii:
        buf_aea = pt_aea.geometry.buffer(radius)
        buf_wgs = gpd.GeoDataFrame(geometry=buf_aea, crs=pt_aea.crs).to_crs("EPSG:4326")
        query_wkt = buf_wgs.geometry.iloc[0].wkt

        q_wkt_mukey = f"""
        SELECT mukey, muname FROM mapunit
        WHERE mukey IN (SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('{query_wkt}'))
        AND mukey IN {_format_in(mukeys_in)}
        """
        mukey_subset = _sda_query(q_wkt_mukey)
        if mukey_subset is not None and not mukey_subset.empty:
            break

    if mukey_subset is None or mukey_subset.empty:
        return str(mukey_df.iloc[0]["mukey"])

    if len(mukey_subset) > 1:
        mukeys_subset_in = mukey_subset["mukey"].dropna().astype(str).unique().tolist()
        q_geom = f"""
        SELECT G.MupolygonWktWgs84 AS geom, mapunit.mukey AS mukey
        FROM mapunit CROSS APPLY SDA_Get_MupolygonWktWgs84_from_Mukey(mapunit.mukey) AS G
        WHERE mukey IN {_format_in(mukeys_subset_in)}
        """
        mukey_geom = _sda_query(q_geom)
        if mukey_geom is not None and not mukey_geom.empty:
            mukey_geom["geometry"] = mukey_geom["geom"].apply(loads)
            geom_gdf = gpd.GeoDataFrame(mukey_geom, geometry="geometry", crs="EPSG:4326")
            geom_gdf = geom_gdf.to_crs(pt_aea.crs)

            # Find closest feature using distances from pt_aea geometry
            distances = geom_gdf.geometry.distance(pt_aea.geometry.iloc[0])
            nearest_idx = distances.idxmin()
            return str(geom_gdf.loc[nearest_idx, "mukey"])

    return str(mukey_subset.iloc[0]["mukey"])


def _aggregate_horizons_to_standard_layers(horizon_tbl: pd.DataFrame, bedrock_depth: float) -> Optional[pd.DataFrame]:
    all_layers = {
        "0-5cm": (0.0, 5.0), "5-20cm": (5.0, 20.0), "20-35cm": (20.0, 35.0),
        "35-50cm": (35.0, 50.0), "50-65cm": (50.0, 65.0), "65-80cm": (65.0, 80.0),
        "80-95cm": (80.0, 95.0), "95-110cm": (95.0, 110.0), "110-125cm": (110.0, 125.0),
        "125-140cm": (125.0, 140.0), "140-155cm": (140.0, 155.0), "155-170cm": (155.0, 170.0),
        "170-185cm": (170.0, 185.0), "185-200cm": (185.0, 200.0)
    }
    valid_layers = {name: limits for name, limits in all_layers.items() if limits[0] < bedrock_depth}
    if valid_layers:
        last_name = list(valid_layers.keys())[-1]
        limits = valid_layers[last_name]
        if limits[1] > bedrock_depth:
            valid_layers[last_name] = (limits[0], bedrock_depth)
    else:
        valid_layers = {f"0-{int(bedrock_depth)}cm": (0.0, bedrock_depth)}

    results_list = []
    cols_to_avg = [
        "dbovendry_r", "dbtenthbar_r", "dbthirdbar_r", "dbfifteenbar_r",
        "wsatiated_r", "wtenthbar_r", "wthirdbar_r", "partdensity", "ksat_r",
        "wfifteenbar_r", "sandtotal_r", "claytotal_r", "silttotal_r", "om_r", "fragvol_r"
    ]

    hz_df = horizon_tbl.copy()
    for col in cols_to_avg + ["hzdept_r", "hzdepb_r"]:
        if col in hz_df.columns:
            hz_df[col] = pd.to_numeric(hz_df[col], errors="coerce")

    for name, (top_d, bot_d) in valid_layers.items():
        # Calculate overlaps
        adj_top = hz_df["hzdept_r"].clip(lower=top_d)
        adj_bottom = hz_df["hzdepb_r"].clip(upper=bot_d)
        thickness = (adj_bottom - adj_top).clip(lower=0)

        overlapping = hz_df[thickness > 0].copy()
        if overlapping.empty:
            continue

        overlapping["thickness"] = thickness[thickness > 0]
        total_thickness = overlapping["thickness"].sum()
        if total_thickness == 0 or np.isnan(total_thickness):
            continue

        summarized = {
            "hzdept_r": top_d,
            "hzdepb_r": bot_d,
            "cokey": overlapping.iloc[0]["cokey"],
            "hzname": "/".join(overlapping["hzname"].dropna().astype(str).unique())
        }

        for col in cols_to_avg:
            if col in overlapping.columns:
                if overlapping[col].isna().all():
                    val = np.nan
                else:
                    val = (overlapping[col] * overlapping["thickness"]).sum() / total_thickness
                summarized[col] = val
            else:
                summarized[col] = np.nan

        results_list.append(summarized)

    if not results_list:
        return None
    return pd.DataFrame(results_list)


def process_soils_ssurgo_alderman(
    grid_points,
    output_dir_csv: str,
    output_dir_individual: str,
    n_cores: int,
    id_col: str,
    lat_col: str,
    long_col: str,
    format_sql_func=None,  # API compatibility; unused
    log_file=None,
    soil_names=None,
    STATSGO: bool = False,
    standardize_layers: bool = False,
) -> bool:
    """
    Query SSURGO using full-profile Alderman logic (dominant component, measured tension
    fallbacks, Saxton & Rawls PTFs) for every point in *grid_points*, write per-point .SOL files.
    """
    print("Starting SSURGO Processing (dominant component/measured tension fallbacks)...")
    _soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "Starting SSURGO Processing (dominant component/measured tension fallbacks)")

    os.makedirs(output_dir_individual, exist_ok=True)

    # Convert GeoDataFrame to df
    gdf = grid_points.copy()
    if hasattr(gdf, "geometry"):
        import geopandas as gpd
        gdf = gdf.to_crs("EPSG:4326")
        gdf[lat_col] = gdf.geometry.y
        gdf[long_col] = gdf.geometry.x

    gdf["_row_idx"] = range(len(gdf))

    if soil_names is None:
        if "soil_name" in gdf.columns:
            soil_names = gdf["soil_name"].tolist()
        elif "soil_names" in gdf.columns:
            soil_names = gdf["soil_names"].tolist()
        elif "compname" in gdf.columns:
            soil_names = gdf["compname"].tolist()

    all_ids = gdf[id_col].astype(str).tolist()

    # Smart resume check
    existing = {
        os.path.splitext(f)[0]
        for f in os.listdir(output_dir_individual)
        if f.endswith(".SOL")
    }

    missing_mask = [str(pid) not in existing for pid in all_ids]
    to_process = gdf[missing_mask].reset_index(drop=True)

    n_total = len(all_ids)
    n_proc = len(to_process)
    n_skip = n_total - n_proc

    print(f"Resume Check: Found {n_skip} existing profiles. Processing {n_proc} remaining.")
    _soil_helper_log(log_file, "INFO", "SSURGO_MAIN", f"Resume check: {n_skip} existing, {n_proc} remaining.")

    if n_proc == 0:
        # Refresh mapping CSV if all exist
        mapping_df = gdf[[id_col, lat_col, long_col]].copy()
        mapping_df = mapping_df.rename(columns={id_col: "ID", lat_col: "latitude", long_col: "longitude"})
        mapping_df["SOIL_ID"] = mapping_df["ID"]
        mapping_df.to_csv(output_dir_csv, index=False)
        print("All soil profiles already exist. Skipping SSURGO processing.")
        _soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "All soil profiles already exist. Skipping SSURGO processing.")
        return True

    CHUNK_SIZE = 1000
    num_chunks = math.ceil(n_proc / CHUNK_SIZE)
    print(f"Processing {n_proc} points in {num_chunks} chunk(s)...")
    _soil_helper_log(log_file, "INFO", "SSURGO_MAIN", f"Processing {n_proc} points in {num_chunks} chunk(s).")

    csv_header_written = os.path.exists(output_dir_csv)

    for chunk_i in range(num_chunks):
        s = chunk_i * CHUNK_SIZE
        e = min((chunk_i + 1) * CHUNK_SIZE, n_proc)
        chunk = to_process.iloc[s:e]
        print(f"  > Chunk {chunk_i+1}/{num_chunks} (Points {s+1} – {e})")
        _soil_helper_log(log_file, "INFO", "SSURGO_CHUNK", f"Chunk {chunk_i+1}/{num_chunks} (Points {s+1}-{e})")

        tasks = [
            {"ID": str(row[id_col]),
             "lat": float(row[lat_col]),
             "lon": float(row[long_col]),
             "output_dir": output_dir_individual,
             "log_file": log_file,
             "soil_name": soil_names[int(row["_row_idx"])] if soil_names is not None else None,
             "STATSGO": STATSGO,
             "standardize_layers": standardize_layers}
            for _, row in chunk.iterrows()
        ]

        results = []
        iter_obj = (
            tqdm(tasks, desc=f"Chunk {chunk_i+1}", unit="pt") if _HAS_TQDM else tasks
        )

        with ThreadPoolExecutor(max_workers=min(n_cores, 16)) as pool:
            future_map = {pool.submit(_process_point, t): t["ID"] for t in tasks}
            for fut in as_completed(future_map):
                pid = future_map[fut]
                try:
                    res = fut.result()
                    if res is not None:
                        results.append(res)
                except Exception as exc:
                    warnings.warn(f"Point {pid} failed: {exc}")
                    _soil_helper_log(log_file, "ERROR", "SSURGO_POINT", f"Processing failed: {exc}", pid)

        if results:
            chunk_df = pd.concat(results, ignore_index=True)

            # Match R column names for final mapping CSV: ID, SOIL_ID, mukey, cokey, compname, comppct_r, latitude, longitude
            # Match the order and merge with base mapping to preserve all IDs
            chunk_df = chunk_df.rename(columns={"ID": id_col})

            # Read or initialize the full mapping CSV
            base_df = gdf[[id_col, lat_col, long_col]].copy()
            base_df = base_df.rename(columns={lat_col: "latitude", long_col: "longitude"})
            base_df["SOIL_ID"] = base_df[id_col].astype(str)

            # Merge new details
            final_mapping = pd.merge(base_df, chunk_df, on=id_col, how="left", suffixes=("", "_new"))
            for col in ["mukey", "cokey", "compname", "comppct_r"]:
                if f"{col}_new" in final_mapping.columns:
                    final_mapping[col] = final_mapping[f"{col}_new"]

            # Keep only required columns
            cols_to_keep = [id_col, "SOIL_ID", "mukey", "cokey", "compname", "comppct_r", "latitude", "longitude"]
            cols_to_keep = [c for c in cols_to_keep if c in final_mapping.columns]
            final_mapping = final_mapping[cols_to_keep].rename(columns={id_col: "ID"})

            final_mapping.to_csv(
                output_dir_csv,
                mode="a" if csv_header_written else "w",
                index=False,
                header=not csv_header_written
            )
            csv_header_written = True

    print("SSURGO Processing Complete.")
    _soil_helper_log(log_file, "INFO", "SSURGO_MAIN", "SSURGO Processing Complete.")
    return True


# ---- Soil Pull by Name / Coordinates API -------------------------------------

def _alderman_coordinates(
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    *,
    long: Optional[float] = None,
    pt_geom=None,
    required: bool = False,
) -> tuple[Optional[float], Optional[float]]:
    """Normalize R/Python coordinate aliases and point-geometry inputs."""
    if lon is not None and long is not None and not math.isclose(float(lon), float(long)):
        raise ValueError("lon and long were both supplied with different values")
    if lon is None:
        lon = long

    if pt_geom is not None:
        geometry = pt_geom
        if hasattr(geometry, "geometry"):
            geometry = geometry.geometry
        if hasattr(geometry, "iloc"):
            if len(geometry) != 1:
                raise ValueError("pt_geom must contain exactly one point")
            geometry = geometry.iloc[0]
        if hasattr(geometry, "geoms") and not hasattr(geometry, "x"):
            geoms = list(geometry.geoms)
            if len(geoms) != 1:
                raise ValueError("pt_geom must contain exactly one point")
            geometry = geoms[0]
        try:
            lon = float(geometry.x)
            lat = float(geometry.y)
        except (AttributeError, TypeError, ValueError) as exc:
            raise TypeError("pt_geom must be a point geometry with x/y coordinates") from exc

    if required and (lat is None or lon is None):
        raise ValueError("Either pt_geom or lat and lon/long must be supplied")
    return (
        None if lat is None else float(lat),
        None if lon is None else float(lon),
    )

def pull_profile_by_name_alderman(
    soil_name: str,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    SSURGO: bool = True,
    STATSGO: bool = False,
    standardize_layers: bool = False,
    log_file: Optional[str] = None,
    profile_id: Optional[str] = None,
    pt_geom=None,
    long: Optional[float] = None,
) -> Optional[dict]:
    """Pull soil profile by name using Alderman SSURGO/STATSGO logic."""
    lat, lon = _alderman_coordinates(lat, lon, long=long, pt_geom=pt_geom)
    _soil_helper_log(log_file, "INFO", "PULL_BY_NAME", f"Pulling profile by name: {soil_name}")

    comp_condition = f"compname IN {_format_in([soil_name])}"

    source_condition = ""
    if not SSURGO:
        source_condition = " AND lgd.areaname = 'United States'"
    elif not STATSGO:
        source_condition = " AND lgd.areaname != 'United States'"

    mukey_query = f"""
    SELECT lgd.areaname AS area_name, compname AS soil_name, co.mukey AS mukey, cokey
    FROM component co
    INNER JOIN mapunit mu ON co.mukey = mu.mukey
    INNER JOIN legend lgd ON mu.lkey = lgd.lkey
    WHERE {comp_condition}{source_condition}
    AND cokey IN (SELECT cokey FROM chorizon)
    """

    mukey_out = _sda_query(mukey_query)
    if mukey_out is None or mukey_out.empty:
        _soil_helper_log(log_file, "WARN", "PULL_BY_NAME", f"No mukeys found for name {soil_name}")
        return None

    if len(mukey_out) > 1 and lat is not None and lon is not None:
        nearest_mukey = _filter_mukey_by_coord_alderman(mukey_out, lat, lon, soil_name)
        mukey_out = mukey_out[mukey_out["mukey"].astype(str) == str(nearest_mukey)]

    comp_cokeys = mukey_out["cokey"].dropna().astype(str).unique().tolist()
    comp_query = f"""
    SELECT compname, cokey, mukey, COALESCE(comppct_r,'') AS comppct_r,
           COALESCE(hydgrp,'') AS hydgrp, COALESCE(slope_r,'') AS slope_r,
           COALESCE(drainagecl,'') AS drainage, COALESCE(albedodry_r,'') AS albedodry_r
    FROM component WHERE cokey IN {_format_in(comp_cokeys)}
    """

    comp_tbl = _sda_query(comp_query)
    if comp_tbl is None or comp_tbl.empty:
        _soil_helper_log(log_file, "WARN", "PULL_BY_NAME", "No component records found for candidate cokeys")
        return None

    comp_row = _choose_dominant_component(comp_tbl)
    horizon_tbl = _get_component_horizons(comp_row.get("cokey"))
    if horizon_tbl is None or horizon_tbl.empty:
        _soil_helper_log(log_file, "WARN", "PULL_BY_NAME", "No horizon records found for dominant component")
        return None

    profile_lat = lat if lat is not None else 39.0
    profile_lon = lon if lon is not None else -95.0
    prof_id = profile_id if profile_id is not None else soil_name

    return _build_dssat_profile_from_component(
        comp_row, horizon_tbl,
        point_id=prof_id,
        lat=profile_lat,
        lon=profile_lon,
        log_file=log_file,
        standardize_layers=standardize_layers
    )


def pull_profile_by_coords_alderman(
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    site: str = "SITE",
    SSURGO: bool = True,
    STATSGO: bool = False,
    standardize_layers: bool = False,
    log_file: Optional[str] = None,
    pt_geom=None,
    long: Optional[float] = None,
) -> Optional[dict]:
    """Pull soil profile by coordinates using Alderman SSURGO/STATSGO logic."""
    lat, lon = _alderman_coordinates(
        lat, lon, long=long, pt_geom=pt_geom, required=True
    )
    _soil_helper_log(log_file, "INFO", "PULL_BY_COORDS", f"Pulling profile by coords: {lat}, {lon}")

    wkt = f"POINT({lon} {lat})"

    source_condition = ""
    if not SSURGO:
        source_condition = " AND lgd.areaname = 'United States'"
    elif not STATSGO:
        source_condition = " AND lgd.areaname != 'United States'"

    mukey_query = f"""
    SELECT lgd.areaname AS area_name, compname AS soil_name, co.mukey AS mukey, cokey
    FROM component co
    INNER JOIN mapunit mu ON co.mukey = mu.mukey
    INNER JOIN legend lgd ON mu.lkey = lgd.lkey
    WHERE mu.mukey IN (
      SELECT * FROM SDA_Get_Mukey_from_intersection_with_WktWgs84('{wkt}')
    ){source_condition}
    AND cokey IN (SELECT cokey FROM chorizon)
    """

    mukey_out = _sda_query(mukey_query)
    if mukey_out is None or mukey_out.empty:
        _soil_helper_log(log_file, "WARN", "PULL_BY_COORDS", "No components found with standard query; trying fallback")
        mukeys = _sda_spatial_mukeys(lat, lon, STATSGO=STATSGO)
        if not mukeys:
            return None
        return _build_simple_fallback_profile(lat, lon, site, mukeys, None, log_file, STATSGO=STATSGO)

    comp_cokeys = mukey_out["cokey"].dropna().astype(str).unique().tolist()
    comp_query = f"""
    SELECT compname, cokey, mukey, COALESCE(comppct_r,'') AS comppct_r,
           COALESCE(hydgrp,'') AS hydgrp, COALESCE(slope_r,'') AS slope_r,
           COALESCE(drainagecl,'') AS drainage, COALESCE(albedodry_r,'') AS albedodry_r
    FROM component WHERE cokey IN {_format_in(comp_cokeys)}
    """

    comp_tbl = _sda_query(comp_query)
    if comp_tbl is None or comp_tbl.empty:
        mukeys = mukey_out["mukey"].dropna().astype(str).unique().tolist()
        return _build_simple_fallback_profile(lat, lon, site, mukeys, None, log_file, STATSGO=STATSGO)

    comp_row = _choose_dominant_component(comp_tbl)
    horizon_tbl = _get_component_horizons(comp_row.get("cokey"))

    profile = None
    if horizon_tbl is not None and not horizon_tbl.empty:
        profile = _build_dssat_profile_from_component(
            comp_row, horizon_tbl,
            point_id=site,
            lat=lat,
            lon=lon,
            log_file=log_file,
            standardize_layers=standardize_layers
        )

    if profile is None:
        mukeys = [str(comp_row.get("mukey"))]
        profile = _build_simple_fallback_profile(lat, lon, site, mukeys, comp_tbl, log_file, STATSGO=STATSGO)

    return profile

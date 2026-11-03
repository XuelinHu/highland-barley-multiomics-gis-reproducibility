#!/usr/bin/env python3
"""Sample all environment rasters at the cleaned GBIF points.

Memory-safe: every raster is sampled point-by-point with rasterio.sample()
(windowed reads under the hood, GDAL block cache).  Points are visited in
latitude order so that striped global rasters (WorldClim/CHELSA) decode each
row strip once.

Outputs
-------
processed/barley_env_matrix.tsv.gz
processed/qingke_plateau_subset.tsv.gz
processed/env_coverage.json
"""
import glob
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import rasterio

PROC = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')
RAS = os.path.join(PROC, "rasters")
RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')

POINTS = os.path.join(PROC, "points_clean.tsv.gz")
SRTM_MOSAIC = os.path.join(RAS, "srtm", "srtm_mosaic.tif")
SLOPE_TIF = os.path.join(RAS, "srtm", "srtm_slope.tif")
ASPECT_TIF = os.path.join(RAS, "srtm", "srtm_aspect.tif")
SOIL_DIR = os.path.join(RAS, "soil")

WC = lambda i: os.path.join(RAW, "worldclim", f"wc2.1_30s_bio_{i}.tif")
WC_ELEV = os.path.join(RAW, "worldclim", "wc2.1_30s_elev.tif")
CH = lambda i: os.path.join(RAW, "chelsa", f"CHELSA_bio{i}_1981-2010_V.2.1.tif")

SOIL_PROPS = ["clay", "sand", "silt", "soc", "phh2o", "nitrogen", "bdod", "cec"]
SOIL_DEPTHS = ["0_5", "5_15"]
# SoilGrids integer -> physical unit scale factors (official SoilGrids convention)
SOIL_SCALE = {
    "clay": 10.0,      # g/kg -> %
    "sand": 10.0,      # g/kg -> %
    "silt": 10.0,      # g/kg -> %
    "soc": 10.0,       # dg/kg -> g/kg
    "phh2o": 10.0,     # pH*10 -> pH
    "nitrogen": 100.0,  # cg/kg -> g/kg
    "bdod": 100.0,     # cg/cm3 -> kg/dm3
    "cec": 10.0,       # mmol(c)/kg -> cmol(c)/kg
}
# CHELSA V2.1 packed-integer decoding (verified against WorldClim v2.1)
CH_KELVIN = [1, 5, 6, 8, 9, 10, 11]   # value/10 - 273.15 -> degC
CH_TEMP = [2, 7]                       # value/10 -> degC
CH_RATIO = [3]                         # value/10 *100 -> %
CH_SEASON = [4]                        # value/10
CH_PREC = [12, 13, 14, 15, 16, 17, 18, 19]  # value/10

LOG = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.append(line)


def sample_sorted(path, lon, lat, order, inv):
    """Sample `path` at (lon, lat); points pre-sorted by `order` (row order)."""
    with rasterio.open(path) as ds:
        coords = zip(lon[order].tolist(), lat[order].tolist())
        vals = np.fromiter((v[0] for v in ds.sample(coords)),
                           dtype=np.float64, count=lon.size)
    return vals[inv]


def sample_soil(lon, lat, order, inv, n, coverage, env):
    """Sample the 16 SoilGrids columns; returns {col: ok}."""
    soil_ok = {}
    for prop in SOIL_PROPS:
        for d in SOIL_DEPTHS:
            # on-disk name keeps the "cm" suffix: sg_clay_0_5cm_mean.tif
            f = os.path.join(SOIL_DIR, f"sg_{prop}_{d}cm_mean.tif")
            col = f"sg_{prop}_{d}"
            env[col] = np.full(n, np.nan)
            if not os.path.exists(f):
                log(f"  !! soil raster missing: {os.path.basename(f)}")
                soil_ok[col] = False
                coverage[col] = {"unit": None, "source": "MISSING RASTER",
                                 "n_valid": 0, "n_nan": int(n), "pct_valid": 0.0,
                                 "min": None, "median": None, "max": None}
                continue
            try:
                v = sample_sorted(f, lon, lat, order, inv)
            except Exception as e:  # noqa: BLE001
                log(f"  !! soil read failed {col}: {type(e).__name__}: {e}")
                soil_ok[col] = False
                coverage[col] = {"unit": None,
                                 "source": f"READ FAILED ({type(e).__name__})",
                                 "n_valid": 0, "n_nan": int(n), "pct_valid": 0.0,
                                 "min": None, "median": None, "max": None}
                continue
            # SoilGrids WCS returns 0 for water / no-data pixels
            v[v <= 0] = np.nan
            v = v / SOIL_SCALE[prop]
            env[col] = v
            soil_ok[col] = True
            cov = {
                "unit": {"clay": "%", "sand": "%", "silt": "%", "soc": "g/kg",
                         "phh2o": "pH", "nitrogen": "g/kg", "bdod": "kg/dm3",
                         "cec": "cmol(c)/kg"}[prop],
                "source": f"SoilGrids 250m {prop} {d.replace('_','-')}cm mean "
                          f"(0 = water/nodata masked)",
                "n_valid": int(np.isfinite(v).sum()),
                "n_nan": int((~np.isfinite(v)).sum()),
                "pct_valid": round(100.0 * np.isfinite(v).sum() / n, 2),
                "min": float(np.nanmin(v)) if np.isfinite(v).any() else None,
                "median": float(np.nanmedian(v)) if np.isfinite(v).any() else None,
                "max": float(np.nanmax(v)) if np.isfinite(v).any() else None,
            }
            coverage[col] = cov
        log(f"  soil {prop} done")
    return soil_ok


def write_outputs(out, coverage, pct_note=""):
    out = out.rename(columns={"lon": "decimalLongitude", "lat": "decimalLatitude"})
    out.to_csv(os.path.join(PROC, "barley_env_matrix.tsv.gz"), sep="\t",
               index=False, compression="gzip")
    sub = out[out.is_qingke_candidate].copy()
    sub.to_csv(os.path.join(PROC, "qingke_plateau_subset.tsv.gz"), sep="\t",
               index=False, compression="gzip")
    log(f"wrote matrix rows={len(out)} cols={out.shape[1]}; "
        f"qingke subset rows={len(sub)} {pct_note}")
    return sub


def soil_only():
    """Re-sample only the 16 soil columns into the existing matrix."""
    t0 = time.time()
    out = pd.read_csv(os.path.join(PROC, "barley_env_matrix.tsv.gz"),
                      sep="\t", low_memory=False)
    n = len(out)
    lon = out.decimalLongitude.to_numpy(dtype=np.float64)
    lat = out.decimalLatitude.to_numpy(dtype=np.float64)
    with rasterio.open(WC(1)) as ds:
        rows = ((ds.bounds.top - lat) / ds.res[0]).astype(np.int64)
    order = np.argsort(rows, kind="stable")
    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)
    log(f"soil-only refresh on {n} points")
    env, coverage = {}, {}
    soil_ok = sample_soil(lon, lat, order, inv, n, coverage, env)
    for k, v in env.items():
        out[k] = v
    sub = write_outputs(out, coverage, "(soil columns refreshed)")
    covj = json.load(open(os.path.join(PROC, "env_coverage.json")))
    covj["coverage"].update(coverage)
    covj["soil_columns_ok"] = soil_ok
    with open(os.path.join(PROC, "env_coverage.json"), "w") as fh:
        json.dump(covj, fh, indent=1, ensure_ascii=False)
    cov = pd.DataFrame(covj["coverage"]).T.reset_index().rename(columns={"index": "variable"})
    cov.to_csv(os.path.join(PROC, "variable_coverage.tsv"), sep="\t", index=False)
    log(f"soil-only done in {time.time()-t0:.0f}s; ok={sum(soil_ok.values())}/16")


def main():
    t0 = time.time()
    pts = pd.read_csv(POINTS, sep="\t", low_memory=False)
    n = len(pts)
    log(f"loaded {n} cleaned points")
    lon = pts.lon.to_numpy(dtype=np.float64)
    lat = pts.lat.to_numpy(dtype=np.float64)

    # visit points in row order (north -> south) for strip-cache friendliness
    with rasterio.open(WC(1)) as ds:
        res = ds.res[0]
        top = ds.bounds.top
    rows = ((top - lat) / res).astype(np.int64)
    order = np.argsort(rows, kind="stable")
    inv = np.empty_like(order)
    inv[order] = np.arange(order.size)

    env = {}
    coverage = {}

    def add(name, arr, unit, source, valid=None):
        env[name] = arr
        v = np.isfinite(arr) if valid is None else valid
        coverage[name] = {
            "unit": unit, "source": source,
            "n_valid": int(v.sum()), "n_nan": int((~v).sum()),
            "pct_valid": round(100.0 * v.sum() / arr.size, 2),
            "min": float(np.nanmin(arr)) if v.any() else None,
            "median": float(np.nanmedian(arr)) if v.any() else None,
            "max": float(np.nanmax(arr)) if v.any() else None,
        }

    # ---------------- terrain -------------------------------------------
    log("sampling SRTM mosaic ...")
    srtm = sample_sorted(SRTM_MOSAIC, lon, lat, order, inv)
    srtm[srtm <= -32000] = np.nan
    log("sampling WorldClim elevation ...")
    welev = sample_sorted(WC_ELEV, lon, lat, order, inv)
    welev[welev <= -32000] = np.nan
    log("sampling slope / aspect ...")
    slope = sample_sorted(SLOPE_TIF, lon, lat, order, inv)
    aspect = sample_sorted(ASPECT_TIF, lon, lat, order, inv)

    elev = np.where(np.isfinite(srtm), srtm, welev)
    elev_source = np.where(np.isfinite(srtm), "SRTM",
                           np.where(np.isfinite(welev), "WorldClim", "none"))
    env["elev_srtm_m"] = srtm
    env["elev_worldclim_m"] = welev
    env["elev"] = elev
    coverage["elev_srtm_m"] = {
        "unit": "m", "source": "SRTM 90m mosaic (CGIAR, 60-110E/20-45N)",
        "n_valid": int(np.isfinite(srtm).sum()), "n_nan": int((~np.isfinite(srtm)).sum()),
        "pct_valid": round(100.0 * np.isfinite(srtm).sum() / n, 2),
        "min": float(np.nanmin(srtm)), "median": float(np.nanmedian(srtm)),
        "max": float(np.nanmax(srtm))}
    coverage["elev_worldclim_m"] = {
        "unit": "m", "source": "WorldClim 2.1 30s elevation",
        "n_valid": int(np.isfinite(welev).sum()), "n_nan": int((~np.isfinite(welev)).sum()),
        "pct_valid": round(100.0 * np.isfinite(welev).sum() / n, 2),
        "min": float(np.nanmin(welev)), "median": float(np.nanmedian(welev)),
        "max": float(np.nanmax(welev))}
    coverage["elev"] = {
        "unit": "m", "source": "SRTM preferred, WorldClim fallback",
        "n_valid": int(np.isfinite(elev).sum()), "n_nan": int((~np.isfinite(elev)).sum()),
        "pct_valid": round(100.0 * np.isfinite(elev).sum() / n, 2),
        "min": float(np.nanmin(elev)), "median": float(np.nanmedian(elev)),
        "max": float(np.nanmax(elev))}
    add("slope", slope, "degrees (0-90)", "SRTM 90m, Horn 3x3")
    add("aspect", aspect, "degrees clockwise from north (0-360), NaN where slope<0.1 deg",
        "SRTM 90m, Horn 3x3")
    log(f"  elev: SRTM {int(np.isfinite(srtm).sum())}, WC "
        f"{int((~np.isfinite(srtm) & np.isfinite(welev)).sum())}, none "
        f"{int((~np.isfinite(elev)).sum())}")

    # ---------------- WorldClim -----------------------------------------
    for i in range(1, 20):
        v = sample_sorted(WC(i), lon, lat, order, inv)
        v[v <= -1e38] = np.nan
        unit = "degC" if i <= 11 else "mm"
        if i == 3:
            unit = "%"
        elif i == 4:
            unit = "index (sd*100)"
        elif i == 15:
            unit = "%"
        add(f"wc_bio{i}", v, unit, f"WorldClim 2.1 30s bio{i} (1970-2000)")
        if i % 5 == 0:
            log(f"  worldclim {i}/19 done")

    wc_land = np.isfinite(env["wc_bio1"])
    log(f"WorldClim land pixels: {int(wc_land.sum())} / {n}")

    # ---------------- CHELSA --------------------------------------------
    for i in range(1, 20):
        v = sample_sorted(CH(i), lon, lat, order, inv)
        v[v >= 65535] = np.nan
        if i in CH_KELVIN:
            v = v / 10.0 - 273.15
            unit = "degC"
        elif i in CH_TEMP:
            v = v / 10.0
            unit = "degC"
        elif i in CH_RATIO:
            v = v / 10.0 * 100.0
            unit = "%"
        elif i in CH_SEASON:
            v = v / 10.0
            unit = "index (sd*100)"
        else:
            v = v / 10.0
            unit = "mm" if i != 15 else "%"
        add(f"ch_bio{i}", v, unit, f"CHELSA V2.1 30s bio{i} (1981-2010)")
        if i % 5 == 0:
            log(f"  chelsa {i}/19 done")

    # ---------------- soil ----------------------------------------------
    soil_ok = sample_soil(lon, lat, order, inv, n, coverage, env)

    # ---------------- assemble ------------------------------------------
    out = pts.copy()
    for k, v in env.items():
        out[k] = v
    out["elev_source"] = elev_source
    out["wc_land"] = wc_land
    out["on_plateau"] = ((out.lon >= 70) & (out.lon <= 105)
                         & (out.lat >= 25) & (out.lat <= 40))
    out["is_qingke_candidate"] = out.on_plateau & (out.elev >= 2000)

    sub = write_outputs(out, coverage)

    summary = {
        "n_points": int(n),
        "n_qingke_candidates": int(len(sub)),
        "n_on_plateau": int(out.on_plateau.sum()),
        "elev_source_counts": {str(k): int(v) for k, v in
                               pd.Series(elev_source).value_counts().items()},
        "n_wc_land": int(wc_land.sum()),
        "soil_columns_ok": soil_ok,
        "coverage": coverage,
    }
    with open(os.path.join(PROC, "env_coverage.json"), "w") as fh:
        json.dump(summary, fh, indent=1, ensure_ascii=False)
    with open(os.path.join(PROC, "extract_env.log"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")

    # compact coverage table
    cov = pd.DataFrame(coverage).T.reset_index().rename(columns={"index": "variable"})
    cov.to_csv(os.path.join(PROC, "variable_coverage.tsv"), sep="\t", index=False)
    log(f"done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    if "--soil-only" in sys.argv:
        sys.exit(soil_only())
    sys.exit(main())

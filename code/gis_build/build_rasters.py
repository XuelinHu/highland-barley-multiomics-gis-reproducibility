#!/usr/bin/env python3
"""Raster preprocessing for the highland barley / GIS environment matrix.

Steps
-----
1. Mosaic the 49 CGIAR SRTM 90 m tiles (60-110E, 20-45N) into one GeoTIFF,
   written tile-by-tile (streaming, never holding the whole mosaic in RAM).
2. Derive slope (degrees) and aspect (degrees clockwise from north) from the
   mosaic with a 3x3 Horn window, in row blocks with a 1-row halo.
3. Mosaic the two SoilGrids longitude blocks (X60-85 + X85-110) per
   property x depth (8 x 2 = 16 outputs).

All outputs are EPSG:4326.  Memory-safe: windowed reads / writes only.
"""
import glob
import json
import os
import sys
import time

import numpy as np
import rasterio
from rasterio.transform import from_origin
from rasterio.windows import Window
from rasterio.merge import merge as rio_merge

RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
PROC = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')
RAS = os.path.join(PROC, "rasters")

SRTM_DIR = os.path.join(RAW, "elevation", "srtm_tif")
SRTM_MOSAIC = os.path.join(RAS, "srtm", "srtm_mosaic.tif")
SLOPE_TIF = os.path.join(RAS, "srtm", "srtm_slope.tif")
ASPECT_TIF = os.path.join(RAS, "srtm", "srtm_aspect.tif")
SOIL_DIR = os.path.join(RAS, "soil")

RES = 1.0 / 1200.0  # 3 arcsec
MOSAIC_BOUNDS = (60.0, 20.0, 110.0, 45.0)  # left, bottom, right, top
SRTM_NODATA = -32768

LOG = []


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    LOG.append(line)


# --------------------------------------------------------------------------
# 1. SRTM mosaic
# --------------------------------------------------------------------------
def build_srtm_mosaic():
    if os.path.exists(SRTM_MOSAIC):
        log(f"SKIP mosaic (exists): {SRTM_MOSAIC}")
        return
    tiles = sorted(glob.glob(os.path.join(SRTM_DIR, "srtm_*.tif")))
    log(f"SRTM mosaic from {len(tiles)} tiles")
    left, bottom, right, top = MOSAIC_BOUNDS
    width = int(round((right - left) / RES))
    height = int(round((top - bottom) / RES))
    transform = from_origin(left, top, RES, RES)
    profile = {
        "driver": "GTiff",
        "dtype": "int16",
        "nodata": SRTM_NODATA,
        "width": width,
        "height": height,
        "count": 1,
        "crs": "EPSG:4326",
        "transform": transform,
        "compress": "LZW",
        "predictor": 2,
        "tiled": True,
        "blockxsize": 512,
        "blockysize": 512,
        "BIGTIFF": "YES",
    }
    log(f"  mosaic grid: {width} x {height} px, res={RES:.10f} deg")
    t0 = time.time()
    with rasterio.open(SRTM_MOSAIC, "w", **profile) as dst:
        for i, tp in enumerate(tiles):
            with rasterio.open(tp) as src:
                b = src.bounds
                col_off = int(round((b.left - left) / RES))
                row_off = int(round((top - b.top) / RES))
                # sanity: tile grid must align with the mosaic grid
                assert abs((b.left - left) / RES - col_off) < 1e-6
                assert abs((top - b.top) / RES - row_off) < 1e-6
                data = src.read(1)
                dst.write(data, 1, window=Window(col_off, row_off, src.width, src.height))
                del data
            if (i + 1) % 10 == 0 or i == len(tiles) - 1:
                log(f"  wrote {i+1}/{len(tiles)} tiles ({time.time()-t0:.0f}s)")
    log(f"  mosaic done in {time.time()-t0:.0f}s -> {os.path.getsize(SRTM_MOSAIC)/1e9:.2f} GB")


# --------------------------------------------------------------------------
# 2. slope / aspect (Horn 3x3)
# --------------------------------------------------------------------------
def build_slope_aspect():
    if os.path.exists(SLOPE_TIF) and os.path.exists(ASPECT_TIF):
        log("SKIP slope/aspect (exist)")
        return
    with rasterio.open(SRTM_MOSAIC) as src:
        W, H = src.width, src.height
        transform = src.transform
        profile = src.profile.copy()
    profile.update(
        dtype="float32", nodata=np.nan, compress="DEFLATE", predictor=3,
        tiled=True, blockxsize=512, blockysize=512, BIGTIFF="YES", count=1,
    )

    # metric spacing
    R = 6371000.0
    deg2m = np.pi / 180.0 * R
    dy = RES * deg2m  # north-south, constant on a sphere
    BLOCK = 1024

    t0 = time.time()
    dst_s = rasterio.open(SLOPE_TIF, "w", **profile)
    dst_a = rasterio.open(ASPECT_TIF, "w", **profile)
    try:
        for r0 in range(0, H, BLOCK):
            r1 = min(H, r0 + BLOCK)
            rstart = max(0, r0 - 1)
            rend = min(H, r1 + 1)
            with rasterio.open(SRTM_MOSAIC) as src:
                z = src.read(1, window=Window(0, rstart, W, rend - rstart)).astype(np.float32)
            if r0 == 0:
                z = np.vstack([z[0:1], z])
            if r1 == H:
                z = np.vstack([z, z[-1:]])
            # edge-pad one column on each side so the 3x3 window output is exactly
            # `width` wide (matches the write window) -- only affects cols 0 / W-1
            z = np.pad(z, ((0, 0), (1, 1)), mode="edge")
            # local index of global row j
            base = r0 - rstart + (1 if r0 == 0 else 0)  # local index of global row r0
            assert base == 1, base

            zs = np.where(z <= -32000, np.nan, z)
            a = zs[0:-2, 0:-2]; b = zs[0:-2, 1:-1]; c = zs[0:-2, 2:]
            d = zs[1:-1, 0:-2]; f = zs[1:-1, 2:]
            g = zs[2:, 0:-2];   h = zs[2:, 1:-1];   i = zs[2:, 2:]

            # east-west spacing varies with latitude (rows of the *output* block)
            gl = np.arange(r0, r1)
            lats = transform.f + (gl + 0.5) * transform.e  # transform.e < 0
            dx = (RES * deg2m * np.cos(np.radians(lats)))[:, None].astype(np.float32)

            dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * dx)
            dzdy_s = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * dy)  # +south

            slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy_s))).astype(np.float32)
            aspect = np.degrees(np.arctan2(-dzdx, dzdy_s))
            aspect = np.mod(aspect, 360.0).astype(np.float32)
            aspect[slope < 0.1] = np.nan  # aspect undefined on (near-)flat ground
            assert slope.shape == (r1 - r0, W), slope.shape
            del a, b, c, d, f, g, h, i, dzdx, dzdy_s
            win = Window(0, r0, W, r1 - r0)
            dst_s.write(slope, 1, window=win)
            dst_a.write(aspect, 1, window=win)
            del z, zs, slope, aspect
            if (r0 // BLOCK) % 5 == 0:
                log(f"  slope/aspect rows {r1}/{H} ({time.time()-t0:.0f}s)")
    finally:
        dst_s.close()
        dst_a.close()
    log(f"  slope/aspect done in {time.time()-t0:.0f}s")


# --------------------------------------------------------------------------
# 3. SoilGrids mosaics (2 longitude blocks -> 1)
# --------------------------------------------------------------------------
def fully_readable(path):
    """True if the raster can be read end-to-end (catches truncated downloads)."""
    try:
        with rasterio.open(path) as ds:
            buf = np.empty((ds.height, ds.width), dtype=ds.dtypes[0])
            ds.read(1, out=buf)
        return True
    except Exception as e:  # noqa: BLE001
        log(f"     read check failed on {os.path.basename(path)}: "
            f"{type(e).__name__}: {e}")
        return False



def build_soil_mosaics():
    props = ["clay", "sand", "silt", "soc", "phh2o", "nitrogen", "bdod", "cec"]
    depths = ["0-5cm", "5-15cm"]
    made = {}
    for p in props:
        for d in depths:
            out = os.path.join(SOIL_DIR, f"sg_{p}_{d.replace('-','_')}_mean.tif")
            made[f"{p}_{d}"] = out
            if os.path.exists(out):
                log(f"SKIP soil {p} {d} (exists)")
                continue
            srcs = [
                os.path.join(RAW, "soilgrids", f"{p}_{d}_mean_X60-85_Y20-45.tif"),
                os.path.join(RAW, "soilgrids", f"{p}_{d}_mean_X85-110_Y20-45.tif"),
            ]
            missing = [s for s in srcs if not os.path.exists(s)]
            if missing:
                log(f"  !! MISSING source for {p} {d}: {missing}")
                continue
            # guard against truncated downloads still being repaired elsewhere
            bad = [os.path.basename(s) for s in srcs if not fully_readable(s)]
            if bad:
                log(f"  !! UNREADABLE (truncated) source for {p} {d}: {bad} -> skipped")
                continue
            t0 = time.time()
            datasets = [rasterio.open(s) for s in srcs]
            try:
                arr, tr = rio_merge(datasets, nodata=-32768, dtype="int16")
            except Exception as e:  # noqa: BLE001
                log(f"  !! merge failed for {p} {d}: {type(e).__name__}: {e}")
                for ds_ in datasets:
                    ds_.close()
                continue
            finally:
                pass
            prof = datasets[0].profile.copy()
            for ds_ in datasets:
                ds_.close()
            prof.update(
                height=arr.shape[1], width=arr.shape[2], transform=tr,
                compress="DEFLATE", predictor=2, tiled=True,
                blockxsize=256, blockysize=256, nodata=-32768, count=1,
            )
            with rasterio.open(out, "w", **prof) as dst:
                dst.write(arr)
            log(f"  soil {p} {d}: {arr.shape[2]}x{arr.shape[1]} in {time.time()-t0:.0f}s "
                f"min={np.nanmin(arr)} max={np.nanmax(arr)} "
                f"zeros={(arr==0).sum()} negfill={(arr<=-30000).sum()}")
            del arr
    return made


def main():
    os.makedirs(os.path.join(RAS, "srtm"), exist_ok=True)
    os.makedirs(SOIL_DIR, exist_ok=True)
    t0 = time.time()
    build_srtm_mosaic()
    build_slope_aspect()
    build_soil_mosaics()
    log(f"ALL RASTER BUILD DONE in {time.time()-t0:.0f}s")
    with open(os.path.join(PROC, "build_rasters.log"), "w") as fh:
        fh.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    sys.exit(main())

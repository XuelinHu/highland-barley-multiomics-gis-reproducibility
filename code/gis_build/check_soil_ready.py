#!/usr/bin/env python3
"""Check that all SoilGrids source blocks are complete; drop derived mosaics
whose sources are (still) truncated so that build_rasters.py rebuilds them."""
import os
import sys

import numpy as np
import rasterio

RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw/soilgrids')
SOIL_DIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed/rasters/soil')
PROPS = ["clay", "sand", "silt", "soc", "phh2o", "nitrogen", "bdod", "cec"]
DEPTHS = ["0-5cm", "5-15cm"]


def readable(path):
    try:
        with rasterio.open(path) as ds:
            buf = np.empty((ds.height, ds.width), dtype=ds.dtypes[0])
            ds.read(1, out=buf)
        return True, ""
    except Exception as e:  # noqa: BLE001
        return False, f"{type(e).__name__}: {e}"


def main():
    bad_src, removed, ready = [], [], []
    for p in PROPS:
        for d in DEPTHS:
            srcs = [os.path.join(RAW, f"{p}_{d}_mean_X60-85_Y20-45.tif"),
                    os.path.join(RAW, f"{p}_{d}_mean_X85-110_Y20-45.tif")]
            problems = []
            for s in srcs:
                if not os.path.exists(s):
                    problems.append(f"{os.path.basename(s)} MISSING")
                    continue
                ok, msg = readable(s)
                if not ok:
                    problems.append(f"{os.path.basename(s)}: {msg}")
            derived = os.path.join(SOIL_DIR, f"sg_{p}_{d.replace('-','_')}_mean.tif")
            if problems:
                bad_src.append((p, d, problems))
                if os.path.exists(derived):
                    os.remove(derived)
                    removed.append(os.path.basename(derived))
                print(f"NOT READY {p} {d}: {'; '.join(problems)}")
            else:
                ok_d = os.path.exists(derived)
                ready.append((p, d, ok_d))
                print(f"ready     {p} {d}  derived={'yes' if ok_d else 'NO (needs build)'}")
    print(f"\nready combos: {len(ready)}/16; unusable sources: {len(bad_src)}")
    if removed:
        print("removed stale derived mosaics:", removed)
    return 1 if bad_src else 0


if __name__ == "__main__":
    sys.exit(main())

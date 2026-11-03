#!/usr/bin/env python3
"""Reconnaissance of all raster inputs: CRS, resolution, size, dtype, nodata, validity."""
import glob
import json
import os
import sys

import rasterio

RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
OUT = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')

groups = {
    "worldclim_bio": sorted(glob.glob(os.path.join(RAW, "worldclim", "wc2.1_30s_bio_*.tif"))),
    "worldclim_elev": sorted(glob.glob(os.path.join(RAW, "worldclim", "wc2.1_30s_elev.tif"))),
    "chelsa": sorted(glob.glob(os.path.join(RAW, "chelsa", "*.tif"))),
    "srtm": sorted(glob.glob(os.path.join(RAW, "elevation", "srtm_tif", "*.tif"))),
    "soilgrids": sorted(glob.glob(os.path.join(RAW, "soilgrids", "*.tif"))),
}

report = {}
for gname, files in groups.items():
    report[gname] = []
    print(f"\n=== {gname}: {len(files)} files ===", flush=True)
    for f in files:
        try:
            with rasterio.open(f) as ds:
                rec = {
                    "file": os.path.basename(f),
                    "crs": str(ds.crs),
                    "epsg": ds.crs.to_epsg() if ds.crs else None,
                    "width": ds.width,
                    "height": ds.height,
                    "res": list(ds.res),
                    "dtype": ds.dtypes[0],
                    "nodata": ds.nodata,
                    "bounds": [round(b, 5) for b in ds.bounds],
                    "count": ds.count,
                    "block": list(ds.block_shapes[0]),
                }
        except Exception as e:  # noqa: BLE001
            rec = {"file": os.path.basename(f), "error": f"{type(e).__name__}: {e}"}
        report[gname].append(rec)

    # print compact summary
    for rec in report[gname][:100]:
        if "error" in rec:
            print(f"  BROKEN {rec['file']}: {rec['error']}", flush=True)
        else:
            print(
                f"  {rec['file']:<45} EPSG={rec['epsg']} {rec['width']}x{rec['height']} "
                f"res={rec['res']} {rec['dtype']} nodata={rec['nodata']} "
                f"bounds={rec['bounds']} block={rec['block']}",
                flush=True,
            )

if not os.path.isdir(OUT):
    os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "recon_rasters.json"), "w") as fh:
    json.dump(report, fh, indent=1)
print("\nWrote recon_rasters.json", flush=True)

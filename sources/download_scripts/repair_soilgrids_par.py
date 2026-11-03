#!/usr/bin/env python3
"""Parallel re-download of truncated SoilGrids coverages with strict validation."""
import subprocess, os, sys, glob, shutil
from concurrent.futures import ThreadPoolExecutor
import rasterio
from rasterio.windows import Window

OUT = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw/soilgrids')
WCS = "https://maps.isric.org/mapserv"
CRS = "http://www.opengis.net/def/crs/EPSG/0/4326"
CHUNKS = {"60-85": ("60", "85"), "85-110": ("85", "110")}

def url_for(prop, depth, x0, x1):
    cov = f"{prop}_{depth}_mean"
    return (f"{WCS}?map=/map/{prop}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage"
            f"&COVERAGEID={cov}&FORMAT=GEOTIFF_INT16&SUBSETTINGCRS={CRS}"
            f"&SUBSET=X({x0},{x1})&SUBSET=Y(20,45)")

def valid(path):
    try:
        with rasterio.open(path) as s:
            w, h = s.width, s.height
            s.read(1, window=Window(w-min(64,w), h-min(256,h), min(64,w), min(256,h)))
            s.read(1, window=Window(0, 0, min(64,w), min(64,h)))
        return True
    except Exception:
        return False

def repair(f):
    base = os.path.basename(f)[:-4]
    prop, depth = base.split("_")[0], base.split("_")[1]
    x0, x1 = CHUNKS[base.split("_X")[1].split("_")[0]]
    for attempt in range(1, 9):
        part = f + ".part"
        if os.path.exists(part): os.remove(part)
        r = subprocess.run(["curl","-s","--max-time","1800","-o",part,url_for(prop,depth,x0,x1)],
                           capture_output=True)
        sz = os.path.getsize(part) if os.path.exists(part) else 0
        if r.returncode == 0 and sz > 1_000_000 and valid(part):
            os.replace(part, f)
            return f"{base}: REPAIRED ({sz})"
        print(f"  {base} attempt {attempt}: curl={r.returncode} size={sz}", flush=True)
    if os.path.exists(f + ".part"): os.remove(f + ".part")
    return f"{base}: FAILED"

bad = [f for f in sorted(glob.glob(f"{OUT}/*.tif")) if not valid(f)]
print(f"broken: {len(bad)} -> {[os.path.basename(b) for b in bad]}", flush=True)
with ThreadPoolExecutor(max_workers=len(bad) or 1) as ex:
    for res in ex.map(repair, bad):
        print(res, flush=True)
print("=== repair pass complete ===", flush=True)

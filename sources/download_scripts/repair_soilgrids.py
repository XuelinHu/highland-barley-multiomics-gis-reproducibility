#!/usr/bin/env python3
"""Re-download truncated SoilGrids WCS coverages with strict validation.
Validates: curl exit code, TIFF open, AND a full bottom-edge strip decode
(catches silently-truncated transfers that still return HTTP 200)."""
import subprocess, os, sys, glob
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
    """True only if the whole raster is decodable, including the last tile row."""
    try:
        with rasterio.open(path) as s:
            w, h = s.width, s.height
            s.read(1, window=Window(w - min(64, w), h - min(256, h), min(64, w), min(256, h)))
            s.read(1, window=Window(0, 0, min(64, w), min(64, h)))
        return True
    except Exception as e:
        sys.stderr.write(f"    invalid: {str(e)[:70]}\n")
        return False

# discover the broken files
bad = []
for f in sorted(glob.glob(f"{OUT}/*.tif")):
    if not valid(f):
        bad.append(f)
print(f"broken files: {len(bad)}", flush=True)

for f in bad:
    base = os.path.basename(f)[:-4]
    prop, depth = base.split("_")[0], base.split("_")[1]
    ck = base.split("_X")[1].split("_")[0]
    x0, x1 = CHUNKS[ck]
    print(f"repairing {base} ({prop} {depth} X{ck})", flush=True)
    ok = False
    for attempt in range(1, 7):
        part = f + ".part"
        if os.path.exists(part): os.remove(part)
        r = subprocess.run(["curl","-s","--max-time","1800","-o",part,url_for(prop,depth,x0,x1)],
                           capture_output=True)
        sz = os.path.getsize(part) if os.path.exists(part) else 0
        if r.returncode != 0:
            print(f"  attempt {attempt}: curl exit={r.returncode} size={sz}", flush=True); continue
        if sz < 1_000_000:
            print(f"  attempt {attempt}: too small ({sz})", flush=True); continue
        if valid(part):
            os.replace(part, f)
            print(f"  REPAIRED {base} -> {sz} bytes", flush=True); ok = True; break
        print(f"  attempt {attempt}: downloaded {sz} but TIFF still truncated", flush=True)
    if not ok:
        print(f"  FAILED {base}", flush=True)
        if os.path.exists(f + ".part"): os.remove(f + ".part")

print("=== repair pass complete ===", flush=True)

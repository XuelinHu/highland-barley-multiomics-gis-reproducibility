#!/usr/bin/env python3
"""Verify downloaded GIS files: size vs server Content-Length + format/CRS integrity."""
import os, sys, subprocess, json, glob, urllib.request

R = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
report = {"ok": [], "warn": [], "fail": []}
def add(kind, label, path, size, msg): report[kind].append({"label":label,"path":path,"size":size,"msg":msg})

def head_len(url):
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent":"Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=60) as r:
            if r.status != 200: return None
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None

def check(url, path, label):
    if not os.path.exists(path):
        add("fail", label, path, 0, "MISSING"); return
    sz = os.path.getsize(path)
    exp = head_len(url) if url else None
    if exp is None:      add("warn", label, path, sz, "no Content-Length from server")
    elif sz == exp:      add("ok",   label, path, sz, f"size matches server ({sz})")
    else:                add("fail", label, path, sz, f"SIZE MISMATCH local={sz} server={exp}")

WC = "https://geodata.ucdavis.edu/climate/worldclim/2_1/base"
check(f"{WC}/wc2.1_30s_bio.zip",  f"{R}/worldclim/wc2.1_30s_bio.zip",  "worldclim/wc2.1_30s_bio.zip")
check(f"{WC}/wc2.1_30s_elev.zip", f"{R}/worldclim/wc2.1_30s_elev.zip", "worldclim/wc2.1_30s_elev.zip")

CB = "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"
for n in range(1, 20):
    f = f"CHELSA_bio{n}_1981-2010_V.2.1.tif"
    check(f"{CB}/{f}", f"{R}/chelsa/{f}", f"chelsa/{f}")

SB = "https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF"
for X in range(49, 59):
    for Y in range(4, 9):
        f = f"srtm_{X:02d}_{Y:02d}.zip"
        check(f"{SB}/{f}", f"{R}/elevation/{f}", f"elevation/{f}")

# zip integrity
for z in sorted(glob.glob(f"{R}/worldclim/*.zip")) + sorted(glob.glob(f"{R}/elevation/*.zip")):
    p = subprocess.run(["unzip","-t",z], capture_output=True, text=True)
    rel = os.path.relpath(z, R); sz = os.path.getsize(z)
    if p.returncode == 0 and "No errors detected" in p.stdout:
        n = p.stdout.count("testing: ")
        add("ok", "zip-integrity", rel, sz, f"unzip -t OK ({n} members)")
    else:
        add("fail", "zip-integrity", rel, sz, "unzip -t FAILED")

# raster integrity via rasterio
try:
    import rasterio
    for t in sorted(glob.glob(f"{R}/chelsa/*.tif")) + sorted(glob.glob(f"{R}/soilgrids/*.tif")):
        rel = os.path.relpath(t, R); sz = os.path.getsize(t)
        try:
            with rasterio.open(t) as s:
                b = s.bounds; w,h,dt,cr = s.width, s.height, s.dtypes[0], str(s.crs)
                s.read(1, out_shape=(1, min(64,h), min(64,w)))   # force a real decode
            add("ok", "raster-integrity", rel, sz,
                f"valid {w}x{h} {dt} {cr} lon[{b.left:.2f},{b.right:.2f}] lat[{b.bottom:.2f},{b.top:.2f}]")
        except Exception as e:
            add("fail", "raster-integrity", rel, sz, f"unreadable: {str(e)[:120]}")
except ImportError:
    add("warn", "raster-integrity", "-", None, "rasterio not available")

# GBIF
g = f"{R}/gbif/gbif_hordeum_vulgare.tsv"
if os.path.exists(g):
    import csv, collections
    lic = collections.Counter(); n = 0; bad = 0
    with open(g, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            n += 1; lic[row.get("license","")] += 1
            try:
                la, lo = float(row["decimalLatitude"]), float(row["decimalLongitude"])
                if not (-90<=la<=90 and -180<=lo<=180): bad += 1
            except Exception: bad += 1
    add("ok", "gbif-rows", "gbif/gbif_hordeum_vulgare.tsv", os.path.getsize(g),
        f"{n} data rows, {bad} invalid coords")
    report["gbif_licenses"] = dict(lic)

print("="*72)
print(f"OK: {len(report['ok'])}   WARN: {len(report['warn'])}   FAIL: {len(report['fail'])}")
print("="*72)
if report["fail"]:
    print("\n### FAILURES"); [print(f"  [{f['label']}] {f['path']}\n      size={f['size']} :: {f['msg']}") for f in report["fail"]]
if report["warn"]:
    print("\n### WARNINGS"); [print(f"  [{w['label']}] {w['path']} :: {w['msg']}") for w in report["warn"]]
print("\n### INTEGRITY DETAIL")
for f in report["ok"]:
    if f["label"] in ("zip-integrity","raster-integrity","gbif-rows"):
        print(f"  [{f['label']}] {f['path']} ({f['size']} B) :: {f['msg']}")
if "gbif_licenses" in report:
    print("\n### GBIF LICENSE BREAKDOWN")
    for k,v in sorted(report["gbif_licenses"].items(), key=lambda x:-x[1]): print(f"  {v:7d}  {k}")
json.dump(report, open(f"{R}/verify_report.json","w"), indent=1, ensure_ascii=False)
print(f"\nJSON -> {R}/verify_report.json")

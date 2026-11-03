#!/usr/bin/env python3
"""Parallel GBIF occurrence download for Hordeum vulgare (taxonKey 2706056) with coordinates."""
import json, time, urllib.request, urllib.parse, os, sys, random, shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

R = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
OUT = f"{R}/gbif/gbif_hordeum_vulgare.tsv"
TMP = f"{R}/gbif/.pages"
TAXON_KEY = 2706056
LIMIT, WORKERS = 300, 6
COLS = ["key","scientificName","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
        "year","month","day","basisOfRecord","occurrenceStatus","countryCode","country",
        "stateProvince","locality","elevation","elevationAccuracy","recordedBy","institutionCode",
        "collectionCode","catalogNumber","datasetKey","publishingCountry","taxonKey","species",
        "issues","establishmentMeans","occurrenceID","gbifID","license"]
BASE = "https://api.gbif.org/v1/occurrence/search"
Q = urllib.parse.urlencode({"taxonKey":TAXON_KEY,"hasCoordinate":"true","hasGeospatialIssue":"false","limit":LIMIT})

def fetch(url, tries=8):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"highland-barley-gis/1.0"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            sys.stderr.write(f"  retry {a+1}/{tries} {type(e).__name__}: {str(e)[:70]}\n"); sys.stderr.flush()
            time.sleep(min(30, 2*(a+1)) + random.random())
    raise RuntimeError("failed: "+url)

def rows_to_file(offset):
    path = f"{TMP}/p{offset:07d}.tsv"
    if os.path.exists(path): return offset, -1
    d = fetch(f"{BASE}?{Q}&offset={offset}")
    res = d.get("results", [])
    with open(path, "w", encoding="utf-8") as fh:
        for r in res:
            vals = []
            for c in COLS:
                v = r.get(c, "")
                if isinstance(v, list): v = ";".join(str(x) for x in v)
                vals.append(str(v).replace("\t"," ").replace("\n"," "))
            fh.write("\t".join(vals)+"\n")
    time.sleep(0.15 + random.random()*0.2)
    return offset, len(res)

first = fetch(f"{BASE}?{Q}&offset=0")
total = first["count"]
print(f"total records: {total}", flush=True)
if os.path.exists(TMP): shutil.rmtree(TMP)
os.makedirs(TMP, exist_ok=True)

# page 0 already fetched in `first`
p0 = f"{TMP}/p{0:07d}.tsv"
with open(p0,"w",encoding="utf-8") as fh:
    for r in first.get("results", []):
        vals=[]
        for c in COLS:
            v=r.get(c,"")
            if isinstance(v,list): v=";".join(str(x) for x in v)
            vals.append(str(v).replace("\t"," ").replace("\n"," "))
        fh.write("\t".join(vals)+"\n")

offsets = list(range(0, total, LIMIT))
print(f"pages: {len(offsets)}  workers: {WORKERS}", flush=True)
done = 0; failed = []
with ThreadPoolExecutor(max_workers=WORKERS) as ex:
    futs = {ex.submit(rows_to_file, o): o for o in offsets}
    for f in as_completed(futs):
        o = futs[f]
        try:
            off, n = f.result(); done += 1
            if done % 40 == 0: print(f"  pages {done}/{len(offsets)}", flush=True)
        except Exception as e:
            failed.append(o); print(f"  PAGE FAILED offset={o}: {str(e)[:80]}", flush=True)

# concatenate in offset order
written = 0
with open(OUT,"w",encoding="utf-8") as out:
    out.write("\t".join(COLS)+"\n")
    for o in offsets:
        p = f"{TMP}/p{o:07d}.tsv"
        if not os.path.exists(p): continue
        with open(p,encoding="utf-8") as fh:
            for line in fh:
                out.write(line); written += 1
print(f"DONE wrote {written} rows -> {OUT}", flush=True)
if failed: print(f"FAILED OFFSETS ({len(failed)}): {failed}", flush=True)
shutil.rmtree(TMP, ignore_errors=True)

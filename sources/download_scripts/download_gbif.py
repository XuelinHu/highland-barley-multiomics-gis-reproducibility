#!/usr/bin/env python3
"""Download GBIF occurrence records for Hordeum vulgare (taxonKey 2706056), with coords."""
import json, time, urllib.request, urllib.parse, os, sys

OUT = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw/gbif/gbif_hordeum_vulgare.tsv')
TAXON_KEY = 2706056
LIMIT = 300
COLS = ["key","scientificName","decimalLatitude","decimalLongitude","coordinateUncertaintyInMeters",
        "year","month","day","basisOfRecord","occurrenceStatus","countryCode","country",
        "stateProvince","locality","elevation","elevationAccuracy","recordedBy","institutionCode",
        "collectionCode","catalogNumber","datasetKey","publishingCountry","taxonKey","species",
        "issues","hasGeospatialIssue","establishmentMeans","occurrenceID","gbifID","license"]

def fetch(url, tries=5):
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent":"highland-barley-gis/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:
            sys.stderr.write(f"  retry {a+1}/{tries} {type(e).__name__}: {e}\n"); sys.stderr.flush()
            time.sleep(2 * (a + 1))
    raise RuntimeError("failed: " + url)

base = "https://api.gbif.org/v1/occurrence/search"
q = urllib.parse.urlencode({"taxonKey":TAXON_KEY,"hasCoordinate":"true","hasGeospatialIssue":"false","limit":LIMIT})
first = fetch(f"{base}?{q}&offset=0")
total = first["count"]
print(f"total records available: {total}", flush=True)

n = 0
with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("\t".join(COLS) + "\n")
    offset = 0
    while offset < total:
        d = first if offset == 0 else fetch(f"{base}?{q}&offset={offset}")
        res = d.get("results", [])
        if not res:
            print(f"  empty page at offset {offset}, stopping", flush=True); break
        for r in res:
            vals = []
            for c in COLS:
                v = r.get(c, "")
                if isinstance(v, list): v = ";".join(str(x) for x in v)
                vals.append(str(v).replace("\t"," ").replace("\n"," "))
            fh.write("\t".join(vals) + "\n")
        n += len(res)
        offset += LIMIT
        if offset % 3000 == 0 or offset >= total:
            print(f"  {n}/{total}", flush=True)
        if d.get("endOfRecords"): break
        time.sleep(0.5)
print(f"DONE wrote {n} rows -> {OUT}", flush=True)

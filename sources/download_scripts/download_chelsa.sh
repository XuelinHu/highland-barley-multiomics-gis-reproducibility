#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# CHELSA v2.1 bioclim 1981-2010, 30 arcsec, global GeoTIFF (new /chelsav2/ path)
OUT=${BARLEY_DATA_ROOT}/gis/raw/chelsa
BASE=https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio
cd "$OUT" || exit 1
for N in $(seq 1 19); do
  f="CHELSA_bio${N}_1981-2010_V.2.1.tif"
  [ -s "$f" ] && { echo "SKIP $f"; continue; }
  wget -c -q --tries=4 --timeout=60 --waitretry=10 -O "$f" "$BASE/$f" \
    && echo "OK $f $(stat -c%s $f)" || echo "FAIL $f"
done
echo "=== CHELSA DONE $(date -Is) ==="

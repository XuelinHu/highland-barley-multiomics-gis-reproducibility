#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# SRTM 4.1 90m region tiles: lon 60-110E (X=49..58), lat 45-20N (Y=4..8)
OUT=${BARLEY_DATA_ROOT}/gis/raw/elevation
BASE=https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF
cd "$OUT" || exit 1
for X in $(seq 49 58); do
  for Y in $(seq 4 8); do
    f=$(printf "srtm_%02d_%02d.zip" $X $Y)
    [ -s "$f" ] && { echo "SKIP $f (exists)"; continue; }
    wget -c -q --tries=4 --timeout=60 --waitretry=5 \
      -O "$f" "$BASE/$f" && echo "OK $f $(stat -c%s $f)" || echo "FAIL $f"
  done
done
echo "=== SRTM DONE ==="

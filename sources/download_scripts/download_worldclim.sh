#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
OUT=${BARLEY_DATA_ROOT}/gis/raw/worldclim
BASE=https://geodata.ucdavis.edu/climate/worldclim/2_1/base
cd "$OUT" || exit 1
for f in wc2.1_30s_bio.zip wc2.1_30s_elev.zip; do
  echo "START $f $(date -Is)"
  wget -c --tries=4 --timeout=60 --waitretry=10 -q -O "$f" "$BASE/$f" \
    && echo "OK $f $(stat -c%s $f)" || echo "FAIL $f"
done
echo "=== WORLDCLIM DONE $(date -Is) ==="

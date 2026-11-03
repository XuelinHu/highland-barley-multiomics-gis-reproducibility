#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# Extract downloaded archives in place (keeps the original .zip files)
R=${BARLEY_DATA_ROOT}/gis/raw
echo "=== unzip WorldClim ==="
for z in "$R"/worldclim/*.zip; do
  [ -s "$z" ] || continue
  echo "-- $(basename "$z")"; unzip -o -q "$z" -d "$R/worldclim/" && echo "   extracted OK"
done
echo "=== unzip SRTM tiles ==="
mkdir -p "$R/elevation/srtm_tif"
n=0
for z in "$R"/elevation/*.zip; do
  [ -s "$z" ] || continue
  unzip -o -q "$z" -d "$R/elevation/srtm_tif/" && n=$((n+1))
done
echo "   extracted $n tiles -> elevation/srtm_tif/"
echo "=== DONE $(date -Is) ==="

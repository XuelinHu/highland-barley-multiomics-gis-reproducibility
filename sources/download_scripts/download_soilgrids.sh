#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# SoilGrids 2.0 via ISRIC WCS, region subset lon 60-110E / lat 20-45N
# Split into 2 longitude chunks to stay under server MAXSIZE=16384
OUT=${BARLEY_DATA_ROOT}/gis/raw/soilgrids
cd "$OUT" || exit 1
WCS="https://maps.isric.org/mapserv"
CRS="http://www.opengis.net/def/crs/EPSG/0/4326"
for prop in clay sand silt soc phh2o nitrogen bdod cec; do
  for depth in 0-5cm 5-15cm; do
    for chunk in "60 85" "85 110"; do
      set -- $chunk; X0=$1; X1=$2
      cov="${prop}_${depth}_mean"
      f="${cov}_X${X0}-${X1}_Y20-45.tif"
      [ -s "$f" ] && { echo "SKIP $f"; continue; }
      for attempt in 1 2 3; do
        curl -s --max-time 900 -o "$f" -w "%{http_code} %{size_download}\n" \
          "$WCS?map=/map/${prop}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=${cov}&FORMAT=GEOTIFF_INT16&SUBSETTINGCRS=${CRS}&SUBSET=X(${X0},${X1})&SUBSET=Y(20,45)" \
          > /tmp/wcs_http_$$.txt 2>&1
        code=$(awk '{print $1}' /tmp/wcs_http_$$.txt); sz=$(awk '{print $2}' /tmp/wcs_http_$$.txt)
        if [ "$code" = "200" ] && [ "$sz" -gt 1000000 ] 2>/dev/null; then
          echo "OK $f $sz"; break
        else
          echo "RETRY($attempt) $f http=$code sz=$sz"; rm -f "$f"; sleep 10
        fi
      done
      [ -s "$f" ] || echo "FAIL $f"
    done
  done
done
echo "=== SOILGRIDS DONE $(date -Is) ==="

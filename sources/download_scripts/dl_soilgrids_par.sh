#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# SoilGrids 2.0 WCS region subsets, parallel (3 concurrent). Writes .part then renames,
# so an existing "$f" always means a COMPLETE file.
OUT=${BARLEY_DATA_ROOT}/gis/raw/soilgrids
cd "$OUT" || exit 1
WCS="https://maps.isric.org/mapserv"
CRS="http://www.opengis.net/def/crs/EPSG/0/4326"
one() {
  prop="$1"; depth="$2"; X0="$3"; X1="$4"
  cov="${prop}_${depth}_mean"; f="${cov}_X${X0}-${X1}_Y20-45.tif"
  [ -s "$f" ] && { echo "SKIP $f"; return 0; }
  for attempt in 1 2 3; do
    rm -f "$f.part"
    code=$(curl -s --max-time 900 -o "$f.part" -w "%{http_code}" \
      "$WCS?map=/map/${prop}.map&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID=${cov}&FORMAT=GEOTIFF_INT16&SUBSETTINGCRS=${CRS}&SUBSET=X(${X0},${X1})&SUBSET=Y(20,45)")
    sz=$(stat -c%s "$f.part" 2>/dev/null || echo 0)
    if [ "$code" = "200" ] && [ "$sz" -gt 1000000 ]; then mv "$f.part" "$f"; echo "OK $f $sz"; return 0; fi
    echo "RETRY($attempt) $f http=$code sz=$sz"; sleep 10
  done
  rm -f "$f.part"; echo "FAIL $f"
}
export -f one; export WCS CRS
for prop in clay sand silt soc phh2o nitrogen bdod cec; do
  for depth in 0-5cm 5-15cm; do
    echo "$prop $depth 60 85"; echo "$prop $depth 85 110"
  done
done | xargs -P 3 -n 4 bash -c 'one "$@"' _
echo "=== SOILGRIDS DONE $(date -Is) ==="

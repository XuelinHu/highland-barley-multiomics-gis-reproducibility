#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
R=${BARLEY_DATA_ROOT}/gis/raw/worldclim
cd "$R" || exit 1
BASE=7285071887
TOTAL=10401038857
echo "truncating main file to $BASE"
truncate -s $BASE wc2.1_30s_bio.zip || exit 1
for p in wc2.1_30s_bio.zip.part0 \
         wc2.1_30s_bio.zip.part1 wc2.1_30s_bio.zip.part1.sub0 wc2.1_30s_bio.zip.part1.sub1 \
         wc2.1_30s_bio.zip.part1.sub2 wc2.1_30s_bio.zip.part1.sub3 \
         wc2.1_30s_bio.zip.part2 wc2.1_30s_bio.zip.part3 wc2.1_30s_bio.zip.part4 wc2.1_30s_bio.zip.part5; do
  [ -f "$p" ] || { echo "MISSING $p"; exit 1; }
  cat "$p" >> wc2.1_30s_bio.zip
  echo "  + $p -> $(stat -c%s wc2.1_30s_bio.zip)"
done
new=$(stat -c%s wc2.1_30s_bio.zip)
if [ "$new" = "$TOTAL" ]; then
  echo "SIZE OK $new"
  rm -f wc2.1_30s_bio.zip.part*
  echo "CLEANED parts"
else
  echo "SIZE FAIL got $new want $TOTAL"
fi
echo "=== ASSEMBLY DONE $(date -Is) ==="

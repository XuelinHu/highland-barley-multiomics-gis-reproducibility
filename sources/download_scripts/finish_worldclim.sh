#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# Finish partially-downloaded ranged parts by re-fetching their REMAINDER with fresh
# parallel connections (the original connections were server-throttled).
OUT=${BARLEY_DATA_ROOT}/gis/raw/worldclim
URL=https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip
cd "$OUT" || exit 1
SUB=4    # sub-connections per remaining range

finish() {   # $1=partfile $2=abs_start $3=abs_end(inclusive)
  local f="$1" s0="$2" e0="$3"
  local have; have=$(stat -c%s "$f" 2>/dev/null || echo 0)
  local start=$(( s0 + have ))
  local want=$(( e0 - s0 + 1 ))
  [ "$have" -ge "$want" ] && { echo "ALREADY $f ($have/$want)"; return 0; }
  local remain=$(( e0 - start + 1 ))
  local chunk=$(( remain / SUB + 1 ))
  echo "FINISH $f: have=$have want=$want refetch=[$start,$e0] ($remain bytes, $SUB conns)"
  local i
  for i in $(seq 0 $((SUB-1))); do
    local a=$(( start + i*chunk )); local b=$(( start + (i+1)*chunk - 1 ))
    [ "$a" -gt "$e0" ] && continue
    [ "$b" -gt "$e0" ] && b="$e0"
    echo "$a $b $f.sub$i"
  done > /tmp/fr_$$.txt
  while read -r a b part; do
    ( curl -s --max-time 3600 -r "$a-$b" -o "$part" "$URL" && echo "  sub $(basename $part) $(stat -c%s "$part")" ) &
  done < /tmp/fr_$$.txt
  wait
  local i2
  for i2 in $(seq 0 $((SUB-1))); do
    [ -f "$f.sub$i2" ] && cat "$f.sub$i2" >> "$f"
  done
  rm -f "$f".sub* /tmp/fr_$$.txt
  local now; now=$(stat -c%s "$f")
  if [ "$now" = "$want" ]; then echo "OK $f $now"; else echo "SHORT $f ($now/$want)"; fi
}

finish wc2.1_30s_bio.zip.part1 7804399716 8323727544
finish wc2.1_30s_bio.zip.part2 8323727545 8843055373

# now assemble the final archive from the 6 ordered parts
total=10401038857
have=$(stat -c%s wc2.1_30s_bio.zip)
echo "assembling: main file currently $have bytes"
i=0
for i in 0 1 2 3 4 5; do
  p="wc2.1_30s_bio.zip.part$i"
  [ -f "$p" ] || { echo "MISSING $p"; continue; }
  cat "$p" >> wc2.1_30s_bio.zip
  echo "  appended $p -> $(stat -c%s wc2.1_30s_bio.zip)"
done
new=$(stat -c%s wc2.1_30s_bio.zip)
if [ "$new" = "$total" ]; then
  echo "OK wc2.1_30s_bio.zip $new (complete)"; rm -f wc2.1_30s_bio.zip.part*
else
  echo "FAIL wc2.1_30s_bio.zip got $new want $total"
fi

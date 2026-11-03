#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# WorldClim downloader with parallel ranged fetch for the remainder of large files.
OUT=${BARLEY_DATA_ROOT}/gis/raw/worldclim
BASE=https://geodata.ucdavis.edu/climate/worldclim/2_1/base
cd "$OUT" || exit 1
PARTS=6

fetch_par() {
  local f="$1"
  local url="$BASE/$f"
  local total
  total=$(curl -sIL --max-time 60 "$url" | grep -i '^content-length:' | tail -1 | tr -d '\r' | cut -d' ' -f2)
  [ -z "$total" ] && { echo "FAIL $f (no content-length)"; return 1; }
  local have
  have=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ "$have" = "$total" ]; then echo "SKIP $f (complete $total)"; return 0; fi
  if [ "$have" -gt "$total" ]; then echo "CORRUPT $f (local $have > remote $total), restarting"; rm -f "$f"; have=0; fi
  echo "RESUME $f from $have / $total"
  local remain=$((total - have)); local chunk=$(( remain / PARTS + 1 ))
  local i=0
  for i in $(seq 0 $((PARTS-1))); do
    local s=$(( have + i*chunk )); local e=$(( have + (i+1)*chunk - 1 ))
    [ "$s" -ge "$total" ] && continue
    [ "$e" -ge "$total" ] && e=$((total-1))
    echo "$s $e $i"
  done > /tmp/ranges_$$.txt
  while read -r s e i; do
    ( curl -s --max-time 3600 -r "$s-$e" -o "$f.part$i" "$url" && \
      echo "part$i done $(stat -c%s "$f.part$i")" ) &
  done < /tmp/ranges_$$.txt
  wait
  # append parts in order
  local ok=1
  for i in $(seq 0 $((PARTS-1))); do
    if [ -f "$f.part$i" ]; then cat "$f.part$i" >> "$f"; fi
  done
  rm -f "$f".part* /tmp/ranges_$$.txt
  local newsz; newsz=$(stat -c%s "$f")
  if [ "$newsz" = "$total" ]; then echo "OK $f $newsz"; else echo "FAIL $f (got $newsz want $total)"; fi
}

for f in wc2.1_30s_bio.zip wc2.1_30s_elev.zip; do fetch_par "$f"; done
echo "=== WORLDCLIM DONE $(date -Is) ==="

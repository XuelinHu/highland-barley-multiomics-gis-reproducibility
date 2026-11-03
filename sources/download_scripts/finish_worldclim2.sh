#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# Finish WorldClim bio with small-chunk parallel fetch + per-chunk retry on fresh connections.
OUT=${BARLEY_DATA_ROOT}/gis/raw/worldclim
URL=https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_bio.zip
CHUNK=20971520      # 20 MB
JOBS=8
WORK=/tmp/wcchunks
cd "$OUT" || exit 1

fetch_one() {   # $1=start $2=end $3=outfile
  local s="$1" e="$2" out="$3"
  local want=$(( e - s + 1 )) a got
  for a in 1 2 3 4 5 6 7 8; do
    curl -s --max-time 240 -r "$s-$e" -o "$out.tmp" "$URL"
    got=$(stat -c%s "$out.tmp" 2>/dev/null || echo 0)
    if [ "$got" = "$want" ]; then mv "$out.tmp" "$out"; return 0; fi
    sleep 2
  done
  rm -f "$out.tmp"; echo "CHUNKFAIL $s $e" >> "$WORK/fails.txt"; return 1
}
export -f fetch_one; export URL WORK

fetch_range() {  # $1=start $2=end $3=outfile
  local start="$1" end="$2" out="$3"
  rm -rf "$WORK"; mkdir -p "$WORK"; : > "$WORK/fails.txt"; : > "$WORK/list.txt"
  local s="$start" n=0 e
  while [ "$s" -le "$end" ]; do
    e=$(( s + CHUNK - 1 )); [ "$e" -gt "$end" ] && e="$end"
    printf "%s %s %s/c%06d\n" "$s" "$e" "$WORK" "$n" >> "$WORK/list.txt"
    s=$(( e + 1 )); n=$(( n + 1 ))
  done
  echo "  fetch_range [$start,$end] -> $out : $n chunks of <=$((CHUNK/1048576))MB, $JOBS jobs"
  xargs -P "$JOBS" -n 3 bash -c 'fetch_one "$@"' _ < "$WORK/list.txt"
  : > "$out"
  local a b f
  while read -r a b f; do cat "$f" >> "$out"; done < "$WORK/list.txt"
  local got; got=$(stat -c%s "$out")
  local want=$(( end - start + 1 ))
  if [ "$got" = "$want" ]; then echo "  OK $out $got"; else echo "  SHORT $out $got/$want"; fi
  [ -s "$WORK/fails.txt" ] && { echo "  chunk failures:"; cat "$WORK/fails.txt"; }
}

echo "== finishing part1 sub0 and sub3 =="
fetch_range 7816933476 7943631993 wc2.1_30s_bio.zip.part1.sub0
fetch_range 8197029030 8323727544 wc2.1_30s_bio.zip.part1.sub3

echo "== finishing part2 remainder =="
fetch_range 8345174201 8843055373 /tmp/part2_rest
if [ -f /tmp/part2_rest ]; then cat /tmp/part2_rest >> wc2.1_30s_bio.zip.part2; rm -f /tmp/part2_rest; fi
echo "  part2 now $(stat -c%s wc2.1_30s_bio.zip.part2) (want 519327829)"

echo "== assembling final archive =="
total=10401038857
have=$(stat -c%s wc2.1_30s_bio.zip)
echo "  main file starts at $have bytes"
for i in 0 1 2 3 4 5; do
  p="wc2.1_30s_bio.zip.part$i"
  [ -f "$p" ] || { echo "  MISSING $p"; continue; }
  cat "$p" >> wc2.1_30s_bio.zip
  echo "  appended $p ($(stat -c%s "$p")) -> $(stat -c%s wc2.1_30s_bio.zip)"
done
new=$(stat -c%s wc2.1_30s_bio.zip)
if [ "$new" = "$total" ]; then
  echo "OK wc2.1_30s_bio.zip $new (COMPLETE)"; rm -f wc2.1_30s_bio.zip.part*
else
  echo "FAIL wc2.1_30s_bio.zip got $new want $total"
fi

echo "== elevation zip (338 MB) with same chunked method =="
EURL=https://geodata.ucdavis.edu/climate/worldclim/2_1/base/wc2.1_30s_elev.zip
ETOTAL=$(curl -sIL --max-time 60 "$EURL" | grep -i '^content-length:' | tail -1 | tr -d '\r' | cut -d' ' -f2)
echo "  elev total=$ETOTAL"
URL="$EURL" fetch_range 0 $((ETOTAL-1)) /tmp/elev_full
if [ -f /tmp/elev_full ] && [ "$(stat -c%s /tmp/elev_full)" = "$ETOTAL" ]; then
  mv /tmp/elev_full wc2.1_30s_elev.zip; echo "OK wc2.1_30s_elev.zip $ETOTAL"
else
  echo "FAIL wc2.1_30s_elev.zip"
fi
echo "=== WORLDCLIM DONE $(date -Is) ==="

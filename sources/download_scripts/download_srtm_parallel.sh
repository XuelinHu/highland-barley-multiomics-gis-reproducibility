#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
# SRTM region tiles, PARALLEL fetch (server throttles per-connection).
# lon 60-110E => X=49..58 ; lat 45-20N => Y=4..8. Tile srtm_49_08 is a known server 404.
OUT=${BARLEY_DATA_ROOT}/gis/raw/elevation
BASE=https://srtm.csi.cgiar.org/wp-content/uploads/files/srtm_5x5/TIFF
cd "$OUT" || exit 1
LIST=$(mktemp)
for X in $(seq 49 58); do for Y in $(seq 4 8); do
  f=$(printf "srtm_%02d_%02d.zip" $X $Y)
  [ "$f" = "srtm_49_08.zip" ] && continue      # known 404
  [ -s "$f" ] && continue                      # already have it
  echo "$f"
done; done > "$LIST"
echo "tiles to fetch: $(wc -l < $LIST)"
fetch_one() {
  f="$1"
  for attempt in 1 2 3; do
    if wget -c -q --tries=2 --timeout=45 --waitretry=3 -O "$f" "$BASE/$f"; then
      echo "OK $f $(stat -c%s $f)"; return 0
    fi
    sleep 5
  done
  echo "FAIL $f"; return 1
}
export -f fetch_one; export BASE
cat "$LIST" | xargs -P 6 -I{} bash -c 'fetch_one "$@"' _ {}
rm -f "$LIST"
echo "=== SRTM DONE $(date -Is) ==="

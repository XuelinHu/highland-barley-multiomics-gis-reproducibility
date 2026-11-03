#!/bin/bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
OUT=${BARLEY_DATA_ROOT}/gis/raw/chelsa
BASE=https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio
cd "$OUT" || exit 1
one() {
  N="$1"; f="CHELSA_bio${N}_1981-2010_V.2.1.tif"
  remote=$(curl -sIL --max-time 60 "$BASE/$f" | grep -i '^content-length:' | tail -1 | tr -d '\r' | cut -d' ' -f2)
  local=$(stat -c%s "$f" 2>/dev/null || echo 0)
  if [ -n "$remote" ] && [ "$local" = "$remote" ]; then echo "SKIP $f ($local)"; return 0; fi
  wget -c -q --tries=4 --timeout=60 --waitretry=5 -O "$f" "$BASE/$f" \
    && echo "OK $f $(stat -c%s "$f")" || echo "FAIL $f"
}
export -f one; export BASE
seq 1 19 | xargs -P 4 -I{} bash -c 'one "$@"' _ {}
echo "=== CHELSA DONE $(date -Is) ==="

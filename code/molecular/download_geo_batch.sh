#!/usr/bin/env bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
set -u
ROOT=${BARLEY_DATA_ROOT}/geo_raw
LOG=$ROOT/download_daemon.log
LOCK=$ROOT/.download.lock
mkdir -p "$ROOT"
exec 9>"$LOCK"
flock -n 9 || exit 0
ACC=(GSE225255 GSE128048 GSE103278 GSE223521 GSE41513 GSE47882 GSE27821 GSE27822 GSE252277)
while true; do
  echo "[$(date '+%F %T')] batch start" >> "$LOG"
  for a in "${ACC[@]}"; do
    d="$ROOT/$a"; mkdir -p "$d"
    p="GSE${a:3:-3}nnn"  # 前 len-3 位数字作桶名：GSE103278→GSE103nnn，GSE27821→GSE27nnn
    for kind in matrix suppl; do
      url="https://ftp.ncbi.nlm.nih.gov/geo/series/$p/$a/$kind/"
      timeout 300 wget -q -r -np -nH --cut-dirs=5 -c --tries=4 --timeout=30 -P "$d" "$url" >> "$LOG" 2>&1 || echo "[$(date '+%F %T')] $a $kind failed/timeout" >> "$LOG"
    done
  done
  echo "[$(date '+%F %T')] batch end; disk usage:" >> "$LOG"
  du -sh "$ROOT"/* >> "$LOG" 2>&1 || true
  sleep 1800
done

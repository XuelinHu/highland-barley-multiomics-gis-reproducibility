#!/usr/bin/env bash
: "${BARLEY_DATA_ROOT:=$PWD/evidence}"
set -euo pipefail
SRC=${BARLEY_DATA_ROOT}/geo_raw
OUT=${BARLEY_DATA_ROOT}/processed
mkdir -p "$OUT/raw_extracted" "$OUT/manifests" "$OUT/metadata" "$OUT/expression"
for tar in "$SRC"/*/*_RAW.tar; do
  [ -f "$tar" ] || continue
  acc=$(basename "$tar" _RAW.tar)
  dest="$OUT/raw_extracted/$acc"
  mkdir -p "$dest"
  tar -xf "$tar" -C "$dest"
  find "$dest" -type f -printf '%f\t%s\n' | sort > "$OUT/manifests/${acc}_files.tsv"
done
find "$SRC" -type f -printf '%p\t%s\n' | sort > "$OUT/manifests/downloaded_files.tsv"
sha256sum "$SRC"/*/*_RAW.tar "$SRC"/*/*.txt.gz "$SRC"/*/*.csv.gz 2>/dev/null > "$OUT/manifests/sha256.txt" || true
echo "Prepared data under $OUT"

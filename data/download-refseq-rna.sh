#!/usr/bin/env bash
# Download all human.*.rna.fna.gz shards from NCBI RefSeq into this directory,
# then gunzip. Skips files that already exist (gz or uncompressed).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

BASE_URL="${BASE_URL:-https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot}"
# Override with: MAX=20 ./download-refseq-rna.sh
MAX="${MAX:-16}"
START="${START:-1}"

echo "Target dir: $SCRIPT_DIR"
echo "Source:     $BASE_URL"
echo "Shards:     $START .. $MAX"
echo

for i in $(seq "$START" "$MAX"); do
  name="human.${i}.rna.fna"
  gz="${name}.gz"
  url="${BASE_URL}/${gz}"

  if [[ -f "$name" ]]; then
    echo "[$i/$MAX] skip (exists): $name"
    continue
  fi

  if [[ -f "$gz" ]]; then
    echo "[$i/$MAX] gunzip existing: $gz"
    gunzip -f "$gz"
    continue
  fi

  echo "[$i/$MAX] download: $gz"
  # -C - resume; -f fail on HTTP errors; -L follow redirects
  if ! curl -fL --retry 3 --retry-delay 2 -C - -o "${gz}.partial" "$url"; then
    echo "  FAILED: $url" >&2
    rm -f "${gz}.partial"
    exit 1
  fi
  mv "${gz}.partial" "$gz"
  echo "[$i/$MAX] gunzip: $gz"
  gunzip -f "$gz"
done

echo
echo "Done. Local *.fna:"
ls -lh human.*.rna.fna 2>/dev/null || ls -lh *.fna 2>/dev/null || true

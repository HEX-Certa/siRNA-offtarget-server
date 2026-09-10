#!/usr/bin/env bash
# Download human.1.rna.fna.gz … human.16.rna.fna.gz from NCBI RefSeq into this
# directory, then gunzip. Skips files that already exist (gz or uncompressed).
# A shard that is 404 on the remote is skipped with a warning (NCBI count varies).
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
  # -C - resume; -L follow redirects. Status via -w (do not use -f: 404 is skip).
  code="$(curl -sSL --retry 3 --retry-delay 2 -C - -o "${gz}.partial" -w "%{http_code}" "$url" || true)"
  if [[ "$code" == "404" ]]; then
    echo "[$i/$MAX] WARN: not on remote, skip: $gz" >&2
    rm -f "${gz}.partial"
    continue
  fi
  if [[ "$code" != "200" && "$code" != "206" ]]; then
    echo "  FAILED: $url (HTTP ${code:-curl})" >&2
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

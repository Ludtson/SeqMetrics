#!/usr/bin/env bash
# ==============================================================================
# Script Name:  obtain_indices.sh
# Description:  Stage 2 Whole Genome Evaluation Engine. Imports computed reference
#               matrices (`cai.coa`, `fop.coa`, `cbi.coa`) and computes genome-wide
#               indices. Zero-Default Failsafe enforced — refuses to fall back to
#               CodonW's built-in generic reference tables.
# Author:       Adekola Owoyemi, Casola Lab, TAMU
# Version:      0.1.0 (pre-release, unshared — see Changelog below for what changed)
# Changelog:    0.1.0 - Removed a duplicated copy of this script's body that had
#                        been accidentally concatenated onto an earlier draft
#                        (two full script bodies, two shebang lines, running
#                        CodonW twice per species with identical arguments).
#                        This is the single, correct version.
# ==============================================================================
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <cleaned_fasta_abs> <species_base_name> <workdir_abs>" >&2
  exit 1
fi

clean_seqfile="$1"
species="$2"
workdir="$3"

if ! command -v codonw >/dev/null 2>&1; then
  echo "ERROR: codonw not found in PATH. Activate the proper conda env." >&2
  exit 1
fi

mkdir -p "${workdir}"
cd "${workdir}"

fasta_local="${species}.fna"
cp "${clean_seqfile}" "${fasta_local}"
log_file="${species}.std"

if [ -f cai.coa ]; then
  echo "[$species] Stage 2: calculating indices with custom references..." >> "${log_file}"
  [ -f fop.coa ] || touch fop.coa
  [ -f cbi.coa ] || touch cbi.coa

  codonw "${fasta_local}" \
    -all_indices \
    -cai_file cai.coa \
    -fop_file fop.coa \
    -cbi_file cbi.coa \
    -nomenu \
    -silent \
    "${species}_indices.out" \
    "${species}_indices.blk" >> "${log_file}" 2>&1
else
  # CRITICAL GATEKEEPER: Script fails hard if default fallback is imminent
  echo "[$species] CRITICAL ERROR: 'cai.coa' not found. Pipeline restricted from using defaults." >> "${log_file}"
  exit 1
fi

if [ ! -s "${species}_indices.out" ]; then
  echo "[$species] ERROR: Stage 2 failed because ${species}_indices.out is empty." >> "${log_file}"
  exit 1
fi

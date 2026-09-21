#!/usr/bin/env bash
# ==============================================================================
# Script Name:  calculate_indices.sh
# Description:  Stage 1 Core Reference Engine. Builds custom CodonW reference
#               matrices (cai.coa, fop.coa, cbi.coa) for one species, using
#               either a biological highly-expressed-gene (HEG) seed set or a
#               statistical top-5%-Fop correspondence-analysis fallback.
# Author:       Adekola Owoyemi, Casola Lab, TAMU
# Version:      0.1.0 (pre-release, unshared — see Changelog below for what changed)
# Changelog:    0.1.0 - Added content validation for built .coa reference files
#                        (line-count sanity check + scan of CodonW's own log
#                        output for its known failure phrases) — existence and
#                        non-empty checks alone do not catch a reference that
#                        CodonW accepted but could not build correctly.
#                      - Added explicit --basis-mode selection (auto/heg/coa).
#                        Previously the choice between HEG-based and COA-based
#                        reference building was implicit: whichever HEG file
#                        happened to exist (or didn't) silently decided it.
#                        That's still the default ("auto") behavior, but it is
#                        now a named, logged choice rather than an accident of
#                        file presence — and "heg" mode hard-fails instead of
#                        silently falling back if no HEG file is found, for
#                        anyone who wants HEG-or-nothing guarantees.
#                      - HEG lookup now checks both naming conventions: the
#                        OrthoFinder-expanded output (`_expanded_heg.txt`) AND
#                        the raw GFF-keyword seed output straight from
#                        extract_heg_ids.py (`_heg_ids.txt`), since orthology
#                        expansion is optional and the two scripts previously
#                        used mismatched filenames with no bridge between them.
# ==============================================================================
set -euo pipefail

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <cleaned_fasta_abs> <species_base_name> <workdir_abs> <heg_dir_or_none> <basis_mode: auto|heg|coa>" >&2
  exit 1
fi

clean_seqfile="$1"
species="$2"
workdir="$3"
heg_dir="$4"
basis_mode="${5:-auto}"

case "${basis_mode}" in
  auto|heg|coa) ;;
  *)
    echo "ERROR: --basis-mode must be one of: auto, heg, coa (got '${basis_mode}')" >&2
    exit 1
    ;;
esac

if ! command -v codonw >/dev/null 2>&1; then
  echo "ERROR: codonw not found in PATH. Activate the proper conda env." >&2
  exit 1
fi

mkdir -p "${workdir}"
cd "${workdir}"

fasta_local="${species}.fna"
cp "${clean_seqfile}" "${fasta_local}"
log_file="${species}.std"

# Check both naming conventions a HEG file might use — orthology-expanded
# first (broader, cross-species-confirmed set), then the raw GFF-keyword
# seeds directly (valid on their own; orthology expansion is optional).
HEG_FILE_EXPANDED="${heg_dir}/${species}_expanded_heg.txt"
HEG_FILE_RAW="${heg_dir}/${species}_heg_ids.txt"

resolve_heg_file() {
  if [[ "${heg_dir}" == "NONE" ]]; then
    return 1
  elif [[ -s "${HEG_FILE_EXPANDED}" ]]; then
    echo "${HEG_FILE_EXPANDED}"
    return 0
  elif [[ -s "${HEG_FILE_RAW}" ]]; then
    echo "${HEG_FILE_RAW}"
    return 0
  else
    return 1
  fi
}

# --- EXPLICIT BASIS-MODE SELECTION ---
case "${basis_mode}" in
  heg)
    if HEG_FILE="$(resolve_heg_file)"; then
      echo "[$species] BASIS MODE (explicit: heg): using ${HEG_FILE}" >> "${log_file}"
      cp "${HEG_FILE}" "high_expression_ids.txt"
    else
      echo "[$species] CRITICAL ERROR: --basis-mode heg was requested but no HEG file found at" >> "${log_file}"
      echo "  ${HEG_FILE_EXPANDED}" >> "${log_file}"
      echo "  or ${HEG_FILE_RAW}" >> "${log_file}"
      echo "Refusing to silently fall back to COA/Fop. Run extract_heg_ids.py (and, optionally," >> "${log_file}"
      echo "map_hegs_by_orthology.py) for this species first, or use --basis-mode auto/coa." >> "${log_file}"
      exit 1
    fi
    ;;
  coa)
    echo "[$species] BASIS MODE (explicit: coa): statistical top-5%-Fop reference, ignoring any HEG file." >> "${log_file}"
    ;;
  auto)
    if HEG_FILE="$(resolve_heg_file)"; then
      echo "[$species] BASIS MODE (auto): found HEG file, using ${HEG_FILE}" >> "${log_file}"
      cp "${HEG_FILE}" "high_expression_ids.txt"
    else
      echo "[$species] BASIS MODE (auto): no HEG file found, defaulting to statistical top-5%-Fop." >> "${log_file}"
    fi
    ;;
esac

# If HEG mode (explicit or auto) didn't already produce the seed list, run
# the statistical COA/Fop fallback to build one.
if [[ ! -s "high_expression_ids.txt" ]]; then
    echo "[$species] Running initial COA to determine Fop..." >> "${log_file}"

    codonw "${fasta_local}" \
      -coa_cu \
      -nomenu \
      -silent \
      "${species}_coa.out" \
      "${species}_coa.blk" >> "${log_file}" 2>&1

    if [ ! -s "${species}_coa.out" ]; then
      echo "[$species] ERROR: Stage 1 failed because ${species}_coa.out is empty." >> "${log_file}"
      exit 1
    fi

    grep -v '^title' "${species}_coa.out" | grep -v '^$' | sort -k8,8nr > "${species}_sorted.out"

    total_genes=$(wc -l < "${species}_sorted.out")
    top_count=$(( (total_genes * 5) / 100 ))
    if [ "${top_count}" -lt 1 ]; then
      top_count=1
    fi

    head -n "${top_count}" "${species}_sorted.out" | awk '{print $1}' > "high_expression_ids.txt"
fi

if [ -s "high_expression_ids.txt" ]; then
  echo "[$species] Extracting sequence fractions..." >> "${log_file}"

  # Direct hash lookup, not regex matching (confirmed ~14 minutes on the real
  # O. sativa run, ~38,489 headers x ~1,924 candidate regexes — an
  # O(genes x candidate-IDs) design). Every ID reaching this point is already
  # a clean "seqN" token produced by remap_seq_ids.py earlier in the
  # pipeline, with no whitespace or special characters to worry about, so an
  # exact hash-lookup does the identical job in O(genes) with no regex at
  # all. This assumption holds specifically because this script is always
  # called on the already-remapped FASTA (see run_species.sh) — if that ever
  # changes, this lookup would need to change with it.
  awk '
    NR==FNR { ids[$1]=1; next }
    /^>/    { matched = (substr($0, 2) in ids) }
    matched { print }
  ' "high_expression_ids.txt" "${fasta_local}" > "high_expression.fna" || true
fi

if [ -s "high_expression.fna" ]; then
  echo "[$species] Stage 1: generating custom reference files..." >> "${log_file}"
  rm -f cai.coa fop.coa cbi.coa

  CODONW_STAGE1_LOG_START=$(wc -l < "${log_file}" 2>/dev/null || echo 0)

  codonw "high_expression.fna" \
    -coa_cu \
    -nomenu \
    -silent \
    "high_subset.out" \
    "high_subset.blk" >> "${log_file}" 2>&1
fi

# --- CONTENT VALIDATION (not just existence) ---
# File-existence + non-empty checks alone don't catch a reference that CodonW
# accepted but couldn't build properly: it can print explicit failure text to
# its own log (e.g. "not calculated", "FAILED", "Problems with the number
# genes") while still emitting a short, technically-non-empty .coa file. Scan
# for those phrases in exactly the log lines CodonW itself just wrote (not
# the whole log, which may contain earlier unrelated content), and refuse to
# treat the reference as valid if any appear — this is a real, confirmed
# CodonW behavior (garbled/incomplete references producing plausible-looking
# but deficient output), not a hypothetical.
if [ -s "high_expression.fna" ]; then
  NEW_LOG_LINES=$(tail -n "+$((CODONW_STAGE1_LOG_START + 1))" "${log_file}")
  if echo "${NEW_LOG_LINES}" | grep -qiE "not calculated|FAILED|Problems with the number genes|no .* found in the high bias dataset"; then
    echo "[$species] CRITICAL ERROR: CodonW reported a reference-building problem" >> "${log_file}"
    echo "  (see the lines above this one in this log). Refusing to trust cai.coa/fop.coa/cbi.coa" >> "${log_file}"
    echo "  even though CodonW may have written non-empty files — a reference built from too" >> "${log_file}"
    echo "  few or unrepresentative genes is not a fallback-safe failure, and CodonW does not" >> "${log_file}"
    echo "  reliably refuse to proceed on its own when this happens." >> "${log_file}"
    rm -f cai.coa fop.coa cbi.coa
  fi
fi

# Minimum-content sanity check: a real cai.coa reference table has one entry
# per sense codon (roughly 60, genetic-code-dependent). A file with only a
# handful of lines exists and is non-empty but is not a usable reference.
MIN_COA_LINES=20
if [ -f cai.coa ]; then
  coa_lines=$(wc -l < cai.coa)
  if [ "${coa_lines}" -lt "${MIN_COA_LINES}" ]; then
    echo "[$species] CRITICAL ERROR: cai.coa has only ${coa_lines} lines (expected >= ${MIN_COA_LINES})." >> "${log_file}"
    echo "  Reference set was likely built from too few genes. Refusing to use it." >> "${log_file}"
    rm -f cai.coa fop.coa cbi.coa
  fi
fi

if [ ! -f cai.coa ]; then
  echo "[$species] WARNING: cai.coa was not produced or was rejected by validation. Stage 2 will halt (no fallback)." >> "${log_file}"
fi

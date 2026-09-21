#!/usr/bin/env bash
# ==============================================================================
# Script Name:  orchestrate_codonw.sh
# Description:  Top-level parallel orchestration engine. Collects inputs, 
#               handles configuration flags (including COA and HEG directories), 
#               manages xargs queues, and consolidates final indices.
# Author:       Adekola Owoyemi, Casola Lab, TAMU
# Version:      0.1.0 (pre-release, unshared — see Changelog below for what changed)
# Changelog:    0.1.0 - Added -b/--basis-mode (auto/heg/coa): makes the choice
#                        between biological (HEG) and statistical (COA/Fop)
#                        reference-building an explicit, logged decision
#                        instead of an implicit consequence of which files
#                        happen to exist. Default 'auto' preserves the
#                        original behavior exactly.
#                      - Added -w/--work-dir: the scratch/sandbox directory
#                        now defaults to <caller's cwd>/work instead of being
#                        anchored to this script's own install location, and
#                        can be overridden explicitly. If a work dir from a
#                        previous run already exists at the resolved path, it
#                        is archived (renamed with a timestamp) rather than
#                        being silently reused or deleted, so prior debugging
#                        artifacts are never lost across separate invocations.
#                        Per-species sandboxes are still hard-reset *within* a
#                        single run (reset_workdir) — that isolation is a
#                        correctness requirement, not just tidiness.
# ==============================================================================
set -euo pipefail

INPUT=""
OUT_BASE_DIR=""
CREATE_MASTER="false"
THREADS=1
COA_DIR="NONE"
HEG_DIR="NONE"
BASIS_MODE="auto"
WORK_DIR_OVERRIDE=""

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SPECIES_SCRIPT="${SCRIPT_DIR}/run_species.sh"

print_usage() {
  cat >&2 <<EOF
Usage:
  $0 -i INPUT -o OUTDIR [-c COA_DIR] [-e HEG_DIR] [-b BASIS_MODE] [-w WORK_DIR] [-m] [-t THREADS]

Options:
  -i, --input       Single FASTA, comma-separated list, or directory.
  -o, --out-dir     Directory for results.
  -c, --coa-dir     Directory containing pre-built .coa profiles (or 'NONE').
  -e, --heg-dir     Directory containing HEG lists — either OrthoFinder-expanded
                    (<species>_expanded_heg.txt) or raw GFF-keyword seeds
                    straight from extract_heg_ids.py (<species>_heg_ids.txt).
                    Both naming conventions are checked. (or 'NONE')
  -b, --basis-mode  How to build the Stage 1 reference set. One of:
                      auto  (default) - use a HEG file if one is found for the
                            species, otherwise fall back to the statistical
                            top-5%-Fop correspondence-analysis method. This is
                            the original implicit behavior, now named and
                            logged explicitly.
                      heg   - require a HEG file for every species; hard-fail
                            (no silent fallback) if one isn't found.
                      coa   - always use the statistical top-5%-Fop method,
                            regardless of whether a HEG file exists.
  -w, --work-dir    Scratch/sandbox directory (default: <cwd>/work). If it
                    already exists from a previous run, it is archived with a
                    timestamp suffix rather than reused or deleted.
  -m, --masterfile  Generate a single master_codonw.tsv.
  -t, --threads     Number of parallel jobs (default: 1).
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)       INPUT="$2"; shift 2 ;;
    -o|--out-dir)     OUT_BASE_DIR="$2"; shift 2 ;;
    -c|--coa-dir)     COA_DIR="$2"; shift 2 ;;
    -e|--heg-dir)     HEG_DIR="$2"; shift 2 ;;
    -b|--basis-mode)  BASIS_MODE="$2"; shift 2 ;;
    -w|--work-dir)    WORK_DIR_OVERRIDE="$2"; shift 2 ;;
    -m|--masterfile)  CREATE_MASTER="true"; shift 1 ;;
    -t|--threads)     THREADS="$2"; shift 2 ;;
    -h|--help)        print_usage; exit 0 ;;
    *)                echo "Unknown option: $1" >&2; print_usage; exit 1 ;;
  esac
done

case "${BASIS_MODE}" in
  auto|heg|coa) ;;
  *) echo "ERROR: --basis-mode must be one of: auto, heg, coa (got '${BASIS_MODE}')" >&2; exit 1 ;;
esac

if [[ -z "${INPUT}" || -z "${OUT_BASE_DIR}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  print_usage
  exit 1
fi

if ! command -v codonw >/dev/null 2>&1; then
  echo "ERROR: codonw not found. Activate your codonw environment." >&2
  exit 1
fi

mkdir -p "${OUT_BASE_DIR}"
OUT_BASE_DIR="$(cd "${OUT_BASE_DIR}" && pwd)"

# Work dir defaults to where the user is running this FROM (their cwd), not
# to this script's own install location. Override with -w/--work-dir.
if [[ -n "${WORK_DIR_OVERRIDE}" ]]; then
  WORK_ROOT="${WORK_DIR_OVERRIDE}"
else
  WORK_ROOT="$(pwd)/work"
fi

# If a work dir from an earlier run already exists at this path, archive it
# (timestamped rename) instead of silently reusing or deleting it — prior
# runs' .std logs and intermediate .coa files are exactly what you need for
# post-mortem debugging, and a plain overwrite destroys them with no warning.
if [[ -d "${WORK_ROOT}" ]]; then
  ARCHIVE_NAME="${WORK_ROOT}_$(date +%Y%m%d_%H%M%S)"
  echo "Existing work dir found at ${WORK_ROOT} — archiving to ${ARCHIVE_NAME}" >&2
  mv "${WORK_ROOT}" "${ARCHIVE_NAME}"
fi
mkdir -p "${WORK_ROOT}"
WORK_ROOT="$(cd "${WORK_ROOT}" && pwd)"

MASTER_FILE="${OUT_BASE_DIR}/master_codonw.tsv"
if [[ "${CREATE_MASTER}" == "true" ]]; then
  rm -f "${MASTER_FILE}"
fi

collect_files() {
  local input="$1"
  local -a paths=()

  if [ -d "${input}" ]; then
    shopt -s nullglob
    for f in "${input}"/*.{fasta,fna,fa,faa}; do
      [ -e "$f" ] || continue
      paths+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
    done
    shopt -u nullglob
  elif [[ "${input}" == *,* ]]; then
    IFS=',' read -ra items <<< "${input}"
    for f in "${items[@]}"; do
      paths+=("$(cd "$(dirname "$f")" && pwd)/$(basename "$f")")
    done
  else
    paths+=("$(cd "$(dirname "${input}")" && pwd)/$(basename "${input}")")
  fi

  for p in "${paths[@]}"; do
    echo "$p"
  done
}

reset_workdir() {
  local dir="$1"

  if [[ -z "${dir}" || "${dir}" == "/" || "${dir}" == "." ]]; then
    echo "ERROR: Refusing to reset unsafe workdir path: '${dir}'" >&2
    exit 1
  fi

  if [[ "${dir}" != "${WORK_ROOT}"/* ]]; then
    echo "ERROR: Refusing to reset path outside WORK_ROOT: '${dir}'" >&2
    exit 1
  fi

  rm -rf "${dir}"
  mkdir -p "${dir}"
}

mapfile -t FILES < <(collect_files "${INPUT}")

if [ "${#FILES[@]}" -eq 0 ]; then
  echo "No valid fasta/fna/fa/faa files found." >&2
  exit 1
fi

echo "Found ${#FILES[@]} file(s) to process."

JOB_LIST_FILE="${WORK_ROOT}/job_queue.$$"
rm -f "${JOB_LIST_FILE}"

# ==============================================================================
# PRE-FLIGHT VALIDATION: Ensure all -c pathways resolve perfectly before queuing
# ==============================================================================
if [[ "${COA_DIR}" != "NONE" ]]; then
  echo "Checking reference profile paths for all target species..."
  missing_profiles=0

  for file in "${FILES[@]}"; do
    base_name="$(basename "${file}")"
    species="${base_name%%.*}"
    
    # Resolve the expected path matching our isolation logic rules
    if [[ "${COA_DIR}" == *"/_coas"* || "${COA_DIR}" == *"_coas/"* || "${COA_DIR}" == *"/${species}"* || "${COA_DIR}" == *"/${species}/"* ]]; then
      CHECK_PATH="${COA_DIR}"
    else
      CHECK_PATH="${COA_DIR}/${species}/${species}_coas"
    fi

    # Verify that the folder exists AND contains the core baseline file
    if [[ ! -d "${CHECK_PATH}" ]] || [[ ! -f "${CHECK_PATH}/cai.coa" ]]; then
      echo "  [CRITICAL MISSING] Species '${species}' requires a valid profile at: ${CHECK_PATH}/cai.coa" >&2
      missing_profiles=$((missing_profiles + 1))
    fi
  done

  if [[ ${missing_profiles} -gt 0 ]]; then
    echo "ERROR: ${missing_profiles} species failed reference validation. Aborting run to prevent silent fallback or worker crashes." >&2
    exit 1
  else
    echo "All reference paths successfully verified. Proceeding to execution."
  fi
fi
# ==============================================================================

# Passes all 6 arguments to the xargs pipeline safely
for file in "${FILES[@]}"; do
  base_name="$(basename "${file}")"
  species="${base_name%%.*}"
  
  # --- TARGETED ISOLATION GATEWAY ---
  if [[ "${COA_DIR}" != "NONE" ]]; then
    # Isolate scoring workspace completely away from previous calibration outputs
    WORKDIR="${WORK_ROOT}/${species}_scoring"
    
    # Check if user pointed to a global base directory or an explicit species folder
    if [[ "${COA_DIR}" == *"/_coas"* || "${COA_DIR}" == *"_coas/"* || "${COA_DIR}" == *"/${species}"* || "${COA_DIR}" == *"/${species}/"* ]]; then
      RESOLVED_COA="${COA_DIR}"
    else
      RESOLVED_COA="${COA_DIR}/${species}/${species}_coas/"
    fi
  else
    # Calibration Mode: Keep original intact pathing completely untouched
    WORKDIR="${WORK_ROOT}/${species}"
    RESOLVED_COA="NONE"
  fi
  # ----------------------------------

  # Enforce clean per-species sandboxes to prevent stale file contamination.
  reset_workdir "${WORKDIR}"
  printf '%s\0%s\0%s\0%s\0%s\0%s\0' "${file}" "${species}" "${WORKDIR}" "${RESOLVED_COA}" "${HEG_DIR}" "${BASIS_MODE}" >> "${JOB_LIST_FILE}"
done

if [[ -s "${JOB_LIST_FILE}" ]]; then
  xargs -0 -n 6 -P "${THREADS}" "${RUN_SPECIES_SCRIPT}" < "${JOB_LIST_FILE}"
else
  echo "No jobs generated. Exiting." >&2
  exit 1
fi

rm -f "${JOB_LIST_FILE}"

# Passes all 6 arguments to the xargs pipeline safely
for file in "${FILES[@]}"; do
  base_name="$(basename "${file}")"
  species="${base_name%%.*}"
  
  # --- MATCH THE WORKDIR TO THE CORRECT RUN MODE ---
  if [[ "${COA_DIR}" != "NONE" ]]; then
    WORKDIR="${WORK_ROOT}/${species}_scoring"
  else
    WORKDIR="${WORK_ROOT}/${species}"
  fi
  # -------------------------------------------------

  SPECIES_OUT="${OUT_BASE_DIR}/${species}"
  mkdir -p "${SPECIES_OUT}"
  
  # Gather metrics safely without copying staging files or polluting outputs
  cp "${WORKDIR}"/*.out "${SPECIES_OUT}/" 2>/dev/null || true
  cp "${WORKDIR}"/*.blk "${SPECIES_OUT}/" 2>/dev/null || true
  cp "${WORKDIR}"/*.std "${SPECIES_OUT}/" 2>/dev/null || true
  
  # Only save raw .coa files if they were newly generated during a Calibration run
  if [[ "${COA_DIR}" == "NONE" ]]; then
    cp "${WORKDIR}"/*.coa "${SPECIES_OUT}/" 2>/dev/null || true
    if [ -d "${WORKDIR}/${species}_coas" ]; then
      cp -r "${WORKDIR}/${species}_coas" "${SPECIES_OUT}/"
    fi
  fi

  if [[ "${CREATE_MASTER}" == "true" ]]; then
    OUT_FILE="${WORKDIR}/${species}_indices.out"
    if [[ -s "${OUT_FILE}" ]]; then
      if [[ ! -f "${MASTER_FILE}" ]]; then
        cat "${OUT_FILE}" > "${MASTER_FILE}"
      else
        grep -v '^title' "${OUT_FILE}" >> "${MASTER_FILE}"
      fi
    fi
  fi
done

if [[ "${CREATE_MASTER}" == "true" ]]; then
  echo "Masterfile written to ${MASTER_FILE}"
fi

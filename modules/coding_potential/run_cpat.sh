#!/usr/bin/env bash
#
# =========================================================================
# CPAT TRAINING AND CLASSIFICATION PIPELINE (run_cpat.sh)
# =========================================================================
#
# DESCRIPTION:
# This script automates the three core phases of the Coding Potential
# Assessment Tool (CPAT) pipeline: Hexamer table generation, Logistic
# Regression Model training, and final transcript classification.
#
# It is designed to ensure a clean workflow by checking dependencies (CPAT 
# executables and Rscript), validating input files, and organizing all 
# intermediate and final results into separate, phase-specific directories.
#
# USAGE:
# ./run_cpat.sh -n <NonCoding_FA> -c <Coding_CDS_FA> -t <Transcripts_FA> -o <PREFIX>
#
# ARGUMENTS:
#   -n: FASTA file containing high-confidence NON-CODING training sequences.
#   -c: FASTA file containing known CODING sequences (CDS) for training.
#   -t: FASTA file containing the TARGET TRANSCRIPTS to be classified.
#   -o: Output PREFIX (e.g., 'Bstricta'). This creates a top-level directory: 
#       <PREFIX>_cpat_output/
#
# DEPENDENCIES:
# - make_hexamer_tab, make_logitModel, cpat (CPAT suite)
# - Rscript (R programming language)
#
# =========================================================================
set -euo pipefail

# --- Configuration & Utility Functions ---

# Function to check if a required input file exists
check_file() {
    if [[ ! -f "$1" ]]; then
        echo "Error: Input file not found: $1" >&2
        exit 1
    fi
}

# Function to check if CPAT executables and R are available in PATH
check_dependencies() {
    local missing=0
    local tools=("make_hexamer_tab" "make_logitModel" "cpat" "Rscript")

    echo "Checking for required executables in PATH..." >&2
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "Error: Dependency '$tool' not found in system PATH. Please install it or update your \$PATH." >&2
            missing=1
        fi
    done

    if [[ "$missing" -eq 1 ]]; then
        exit 1
    fi
    echo "All CPAT tools and Rscript found." >&2
}

# --- Argument Parsing ---

while getopts "n:c:t:o:h" opt; do
  case "$opt" in
    n) NONCODING_FA="$OPTARG" ;;     # non-coding training FASTA
    c) CODING_CDS_FA="$OPTARG" ;;    # coding CDS training FASTA
    t) TRANSCRIPTS_FA="$OPTARG" ;;   # transcripts to classify
    o) PREFIX="$OPTARG" ;;           # output prefix (species tag)
    h)
      echo "Usage: $0 -n noncoding_train.fa -c coding_cds.fa -t transcripts.fa -o PREFIX"
      echo ""
      echo "Outputs are organized under: PREFIX_cpat_output/"
      exit 0
      ;;
    *)
      echo "Usage: $0 -n noncoding_train.fa -c coding_cds.fa -t transcripts.fa -o PREFIX"
      exit 1
      ;;
  esac
done

# Check required args
if [[ -z "${NONCODING_FA:-}" || -z "${CODING_CDS_FA:-}" || -z "${TRANSCRIPTS_FA:-}" || -z "${PREFIX:-}" ]]; then
  echo "Error: missing required arguments."
  echo "Usage: $0 -n noncoding_train.fa -c coding_cds.fa -t transcripts.fa -o PREFIX"
  exit 1
fi

# --- Pre-Execution Checks ---
check_dependencies # Checks for CPAT and Rscript
check_file "${NONCODING_FA}"
check_file "${CODING_CDS_FA}"
check_file "${TRANSCRIPTS_FA}"

# --- Setup Output Directory Structure ---
OUTPUT_ROOT="${PREFIX}_cpat_output"
HEXAMER_DIR="${OUTPUT_ROOT}/1_hexamer"
MODEL_DIR="${OUTPUT_ROOT}/2_model"
CLASSIFY_DIR="${OUTPUT_ROOT}/3_classification"

echo "Creating output directories under: ${OUTPUT_ROOT}" >&2
mkdir -p "${HEXAMER_DIR}" "${MODEL_DIR}" "${CLASSIFY_DIR}"

# --- Define Output Paths (using the new directories) ---
HEXAMER_TSV="${HEXAMER_DIR}/${PREFIX}_Hexamer.tsv"
LOG_STEP1="${HEXAMER_DIR}/${PREFIX}_step1.log"

LOGIT_MODEL_PATH="${MODEL_DIR}/${PREFIX}" # Path for model output prefix
LOGIT_MODEL="${LOGIT_MODEL_PATH}.logit.RData"
LOG_STEP2="${MODEL_DIR}/${PREFIX}_step2.log"

CPAT_OUTPUT_PREFIX="${CLASSIFY_DIR}/${PREFIX}_cpat_results"
LOG_STEP3="${CLASSIFY_DIR}/${PREFIX}_step3.log"

# ---

## 💻 PHASE 1: Build Hexamer Table (`make_hexamer_tab`)

echo "=== Phase 1: Build Hexamer Table ===" >&2
make_hexamer_tab \
  -c "${CODING_CDS_FA}" \
  -n "${NONCODING_FA}" \
  > "${HEXAMER_TSV}" 2> "${LOG_STEP1}"

echo "Hexamer table generated: ${HEXAMER_TSV}" >&2

# ---

## 📈 PHASE 2: Build Species-Specific Logit Model (`make_logitModel`)

echo "=== Phase 2: Build Logit Model ===" >&2
# NOTE: The -o argument here accepts a PREFIX for the feature and RData files.
make_logitModel \
  -x "${HEXAMER_TSV}" \
  -c "${CODING_CDS_FA}" \
  -n "${NONCODING_FA}" \
  -o "${LOGIT_MODEL_PATH}" 2> "${LOG_STEP2}"

echo "Logit model trained: ${LOGIT_MODEL}" >&2

# ---

## 📊 PHASE 3: Run CPAT on Target Transcripts (`cpat`)

echo "=== Phase 3: Run CPAT Classification ===" >&2
# NOTE: The -o argument here accepts an OUTPUT PREFIX for the final result files.
cpat \
  -g "${TRANSCRIPTS_FA}" \
  -d "${LOGIT_MODEL}" \
  -x "${HEXAMER_TSV}" \
  -o "${CPAT_OUTPUT_PREFIX}" 2> "${LOG_STEP3}"

echo "Classification complete. Results are in ${CLASSIFY_DIR}/" >&2
echo "Done. All results are organized under the ${OUTPUT_ROOT} directory." >&2

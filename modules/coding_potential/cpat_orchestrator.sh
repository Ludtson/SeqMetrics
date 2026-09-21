#!/usr/bin/env bash
set -euo pipefail

# cpat_orchestrator.sh
# Orchestrates noncoding generation + CPAT training/classification.
#
# Modes:
#   1) Single species (each flag is a file):
#      --genome g.fa --protein p.fa --gff a.gff3 --cds c.fa --out-dir OUT
#
#   2) Multiple species, comma-separated lists (matching by index):
#      --genome "g1.fa,g2.fa" --protein "p1.fa,p2.fa" \
#      --gff "a1.gff3,a2.gff3" --cds "c1.fa,c2.fa" --out-dir OUT
#
#   3) Directory mode (CDS-driven; other dirs must contain matching basenames):
#      --genome /dir/genomes --protein /dir/proteins \
#      --gff /dir/gffs --cds /dir/cds --out-dir OUT
#
# Requirements:
#   - get_noncoding_data.py on PATH
#   - run_cpat.sh on PATH
#   - CPAT + Python available in your env

GENOME_IN=""
PROTEIN_IN=""
GFF_IN=""
CDS_IN=""
OUT_DIR=""

print_usage() {
  cat >&2 <<EOF
Usage:
  $0 --genome GENOME.fa[,...|DIR] --protein PROTEIN.fa[,...|DIR] \\
     --gff ANNOT.gff3[,...|DIR] --cds CDS.fa[,...|DIR] --out-dir OUTDIR

Notes:
- In list mode, comma-separated lists must have the same length.
- In directory mode, CDS dir (--cds) drives the loop, and other dirs
  must contain files whose basenames match the CDS basenames.
EOF
}

# --------- parse args ---------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --genome)  GENOME_IN="$2";  shift 2 ;;
    --protein) PROTEIN_IN="$2"; shift 2 ;;
    --gff)     GFF_IN="$2";     shift 2 ;;
    --cds)     CDS_IN="$2";     shift 2 ;;
    --out-dir) OUT_DIR="$2";    shift 2 ;;
    -h|--help)
      print_usage
      exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      print_usage
      exit 1 ;;
  esac
done

if [[ -z "${GENOME_IN}" || -z "${PROTEIN_IN}" || -z "${GFF_IN}" || -z "${CDS_IN}" || -z "${OUT_DIR}" ]]; then
  echo "ERROR: Missing required arguments." >&2
  print_usage
  exit 1
fi

# --------- dependency checks ---------
if ! command -v get_noncoding_data.py >/dev/null 2>&1; then
  echo "ERROR: get_noncoding_data.py not found in PATH." >&2
  exit 1
fi

if ! command -v run_cpat.sh >/dev/null 2>&1; then
  echo "ERROR: run_cpat.sh not found in PATH." >&2
  exit 1
fi

mkdir -p "${OUT_DIR}"

WKDIR_ROOT="${HOME}/wkdir_cpat"
mkdir -p "${WKDIR_ROOT}"

# --------- helpers ---------
expand_list() {
  local spec="$1"
  local -n out_arr=$2
  IFS=',' read -ra out_arr <<< "$spec"
}

# Find file in DIR whose basename stem matches SPECIES (stem = name before first dot)
find_matching_by_stem() {
  local dir="$1"
  local species="$2"
  local ext_pattern="${3:-*}"  # optional, e.g., "*.gff3" or "*.fa"
  local stem
  local f

  if [[ ! -d "$dir" ]]; then
    echo "ERROR: Expected directory: $dir" >&2
    return 1
  fi

  shopt -s nullglob
  for f in "$dir"/*; do
    [[ -f "$f" ]] || continue
    local bn
    bn=$(basename "$f")
    stem="${bn%%.*}"
    if [[ "$stem" == "$species" ]]; then
      echo "$f"
      shopt -u nullglob
      return 0
    fi
  done
  shopt -u nullglob
  return 1
}

# --------- decide mode ---------
if [[ -d "${CDS_IN}" ]]; then
  MODE="dir"
else
  MODE="list"
fi

echo "Running in MODE=${MODE}"

# --------- directory mode (CDS-driven) ---------
if [[ "${MODE}" == "dir" ]]; then
  CDS_DIR="${CDS_IN}"

  if [[ ! -d "${GENOME_IN}" || ! -d "${PROTEIN_IN}" || ! -d "${GFF_IN}" ]]; then
    echo "ERROR: In directory mode, --genome, --protein, --gff, and --cds must all be directories." >&2
    exit 1
  fi

  shopt -s nullglob
  CDS_FILES=("${CDS_DIR}"/*)
  shopt -u nullglob

  if [[ ${#CDS_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No files found in CDS directory: ${CDS_DIR}" >&2
    exit 1
  fi

  echo "Found ${#CDS_FILES[@]} CDS files to process (CDS-driven)."

  for cds_path in "${CDS_FILES[@]}"; do
    cds_abs=$(realpath "${cds_path}")
    cds_bn=$(basename "${cds_abs}")
    species="${cds_bn%%.*}"

    genome_match=$(find_matching_by_stem "${GENOME_IN}" "${species}") || {
      echo "ERROR: No genome match for species '${species}' in ${GENOME_IN}" >&2
      exit 1
    }

    protein_match=$(find_matching_by_stem "${PROTEIN_IN}" "${species}") || {
      echo "ERROR: No protein match for species '${species}' in ${PROTEIN_IN}" >&2
      exit 1
    }

    gff_match=$(find_matching_by_stem "${GFF_IN}" "${species}") || {
      echo "ERROR: No GFF match for species '${species}' in ${GFF_IN}" >&2
      exit 1
    }

    genome_abs=$(realpath "${genome_match}")
    protein_abs=$(realpath "${protein_match}")
    gff_abs=$(realpath "${gff_match}")

    echo "------------------------------------------------------"
    echo "Species: ${species}"
    echo "  Genome : ${genome_abs}"
    echo "  Protein: ${protein_abs}"
    echo "  GFF    : ${gff_abs}"
    echo "  CDS    : ${cds_abs}"

    WORKDIR="${WKDIR_ROOT}/${species}"
    mkdir -p "${WORKDIR}"
    SPECIES_OUT="${OUT_DIR}/cpat/${species}"
    mkdir -p "${SPECIES_OUT}"

    echo "[${species}] Copying inputs into workdir: ${WORKDIR}"
    cp "${genome_abs}"  "${WORKDIR}/genome.fa"
    cp "${protein_abs}" "${WORKDIR}/proteins.fa"
    cp "${gff_abs}"     "${WORKDIR}/annot.gff3"
    cp "${cds_abs}"     "${WORKDIR}/cds.fa"

    cd "${WORKDIR}"

    echo "[${species}] Generating noncoding data..."

    get_noncoding_data.py \
      -p proteins.fa \
      -g genome.fa \
      -a annot.gff3 \
      -o nc_out \
      -l 200 \
      -m 250 \
      -r "${species}" \
      -n 0

    protein_base=$(basename "proteins.fa")
    protein_stem="${protein_base%.*}"
    protein_ext=".${protein_base##*.}"

    NONCODING_FA="nc_out/${protein_stem}_nc${protein_ext}"

    echo "[${species}] Noncoding FASTA expected at: ${NONCODING_FA}"

    if [[ ! -f "${NONCODING_FA}" ]]; then
      echo "ERROR: Noncoding output not found: ${NONCODING_FA}" >&2
      exit 1
    fi

    echo "[${species}] Running CPAT with run_cpat.sh"

    run_cpat.sh \
      -n "${NONCODING_FA}" \
      -c "cds.fa" \
      -t "cds.fa" \
      -o "${species}"

    echo "[${species}] Copying outputs from workdir to ${SPECIES_OUT}"
    cp -r "${WORKDIR}"/* "${SPECIES_OUT}/"

    echo "[${species}] Done."
  done

  echo "All species completed. Results under: ${OUT_DIR}/cpat/"
  exit 0
fi

# --------- list/single-file mode (index-based) ---------
declare -a GENOMES_ARR PROTEINS_ARR GFFS_ARR CDS_ARR

expand_list "${GENOME_IN}"  GENOMES_ARR
expand_list "${PROTEIN_IN}" PROTEINS_ARR
expand_list "${GFF_IN}"     GFFS_ARR
expand_list "${CDS_IN}"     CDS_ARR

n_g=${#GENOMES_ARR[@]}
n_p=${#PROTEINS_ARR[@]}
n_a=${#GFFS_ARR[@]}
n_c=${#CDS_ARR[@]}

if [[ $n_g -ne $n_p || $n_g -ne $n_a || $n_g -ne $n_c ]]; then
  echo "ERROR: In list mode, --genome, --protein, --gff, and --cds lists must have the same length." >&2
  exit 1
fi

echo "Found ${n_g} species to process (list/single-file mode)."

for idx in "${!GENOMES_ARR[@]}"; do
  genome_abs=$(realpath "${GENOMES_ARR[$idx]}")
  protein_abs=$(realpath "${PROTEINS_ARR[$idx]}")
  gff_abs=$(realpath "${GFFS_ARR[$idx]}")
  cds_abs=$(realpath "${CDS_ARR[$idx]}")

  cds_bn=$(basename "${cds_abs}")
  species="${cds_bn%%.*}"

  echo "------------------------------------------------------"
  echo "Species index ${idx}: ${species}"
  echo "  Genome : ${genome_abs}"
  echo "  Protein: ${protein_abs}"
  echo "  GFF    : ${gff_abs}"
  echo "  CDS    : ${cds_abs}"

  WORKDIR="${WKDIR_ROOT}/${species}"
  mkdir -p "${WORKDIR}"
  SPECIES_OUT="${OUT_DIR}/cpat/${species}"
  mkdir -p "${SPECIES_OUT}"

  echo "[${species}] Copying inputs into workdir: ${WORKDIR}"
  cp "${genome_abs}"  "${WORKDIR}/genome.fa"
  cp "${protein_abs}" "${WORKDIR}/proteins.fa"
  cp "${gff_abs}"     "${WORKDIR}/annot.gff3"
  cp "${cds_abs}"     "${WORKDIR}/cds.fa"

  cd "${WORKDIR}"

  echo "[${species}] Generating noncoding data..."

  get_noncoding_data.py \
    -p proteins.fa \
    -g genome.fa \
    -a annot.gff3 \
    -o nc_out \
    -r "${species}" 

    
  protein_base=$(basename "proteins.fa")
  protein_stem="${protein_base%.*}"
  protein_ext=".${protein_base##*.}"

  NONCODING_FA="nc_out/${protein_stem}_nc${protein_ext}"

  echo "[${species}] Noncoding FASTA expected at: ${NONCODING_FA}"

  if [[ ! -f "${NONCODING_FA}" ]]; then
    echo "ERROR: Noncoding output not found: ${NONCODING_FA}" >&2
    exit 1
  fi

  echo "[${species}] Running CPAT with run_cpat.sh"

  run_cpat.sh \
    -n "${NONCODING_FA}" \
    -c "cds.fa" \
    -t "cds.fa" \
    -o "${species}"
  
  BEST_ORF="${species}cpatoutput/3classification/${species}cpatresults.ORF_prob.best.tsv"
  REPAIRED_BEST_ORF="${species}cpatoutput/3classification/${species}cpatresults.ORF_prob.best.repaired.tsv"

  repair_cpat_ids.py \
    --fasta cds.fa \
    --cpat-table "$BEST_ORF" \
    --out "$REPAIRED_BEST_ORF"  

  echo "[${species}] Cleaning workdir (removing large input files)..."
  rm -f genome.fa cds.fa proteins.fa annot.gff3

  echo "[${species}] Copying outputs from workdir to ${SPECIES_OUT}"
  cp -r "${WORKDIR}"/* "${SPECIES_OUT}/"

  echo "[${species}] Done."
done

echo "All species completed. Results under: ${OUT_DIR}/cpat/"

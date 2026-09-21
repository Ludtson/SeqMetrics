# pepstats-pipeline (v1.0.0)

## Overview

pepstats-pipeline is a minimal, portable workflow for running EMBOSS `pepstats`
on protein FASTA files and converting the results into structured TSV tables
for downstream analysis.

The design prioritizes:

- Simplicity
- Portability (no internal PATH assumptions)
- Reproducibility
- Minimal dependencies

---

## Project Structure

```
pepstats-pipeline/
├── bin/
│   ├── pepstats-run      # Bash pipeline (main entry point)
│   ├── pepstats-parse    # Python parser (standard library only)
│   └── pepflow           # Optional CLI wrapper
│
├── env/
│   └── environment.yml
│
├── results/              # Created during execution
│   ├── pepstats_raw/
│   ├── tables/
│   ├── master/
│   └── logs/
│
└── README.md
```

---

## Requirements

- Conda (Miniconda or Anaconda)
- EMBOSS (`pepstats`)
- Python 3

---

## Environment Setup

From the project root directory:

```bash
cd pepstats-pipeline

conda env create -f env/environment.yml
conda activate pepstats_env
```

Verify that `pepstats` is available:

```bash
which pepstats
```

If not installed:

```bash
conda install -c bioconda emboss
```

---

## EMBOSS Configuration (Important)

EMBOSS tools (including `pepstats`) rely on internal configuration (ACD files)
and data files. In some Conda environments, these are not automatically detected.

If you encounter errors such as:

```
ACD file not opened
```

or:

```
Unable to open data file 'Eamino.dat'
```

you need to explicitly define the required environment variables.

---

### Permanent Fix (Recommended)

Run the following commands once:

```bash
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
nano $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```

Add the following lines:

```bash
export EMBOSS_ACDROOT=$CONDA_PREFIX/share/EMBOSS/acd
export EMBOSS_DATA=$CONDA_PREFIX/share/EMBOSS/data
```

Save the file, then reload the environment:

```bash
conda deactivate
conda activate pepstats-pipeline
```

---

### Verify Configuration

```bash
echo $EMBOSS_ACDROOT
echo $EMBOSS_DATA
```

Then test:

```bash
pepstats -help
```

If the help message prints without errors, the configuration is correct.

---

### Notes

- `EMBOSS_ACDROOT` points to tool definition files (`*.acd`)
- `EMBOSS_DATA` provides required runtime data (e.g. `Eamino.dat`)
- These variables must be set for EMBOSS tools to function correctly in Conda environments

Save the file, then reload the environment:

```bash
conda deactivate
conda activate pepstats-pipeline

```

### Verify

```bash
pepstats -help
```

If the help message prints, the configuration is correct.

---

## Make Scripts Executable

```bash
chmod +x bin/*
```

---

## Running the Pipeline

Run all commands from the project root.

### Single FASTA file

```bash
./bin/pepstats-run -f data/proteins.fa -o output/
```

### Directory of FASTA files

```bash
./bin/pepstats-run -f data/ -o output/
```

### Multiple files (comma-separated)

```bash
./bin/pepstats-run -f a.fa,b.fa -o output/
```

### Generate master table

```bash
./bin/pepstats-run -f data/ -o output/ -m
```

---

## Optional Wrapper (pepflow)

The wrapper provides a simple CLI with subcommands:

```bash
./bin/pepflow run -f data/ -o output/
./bin/pepflow parse -i file.pepstats -o file.tsv
```

---

## Command-Line Options

| Flag | Description |
|------|-------------|
| -f, --fasta   | Input FASTA file, directory, or comma-separated list |
| -o, --outdir  | Output directory |
| -m, --master  | Generate combined master TSV file |
| -d, --dayhoff | Include Dayhoff statistics |
| -t, --trimmed | Output compact feature table |
| -h, --help    | Show help message |

---

## Output Structure

All results are written inside the directory provided with `-o`.

```
output/
├── pepstats_raw/
├── tables/
├── master/
└── logs/
```

---

## Output Format

Each TSV file contains protein statistics and derived features including:

- Gene_ID (sequence identifier)
- Residues (sequence length)
- Average_Residue_Weight
- Charge
- Isoelectric_Point
- Charge_per_residue
- Basic_minus_Acidic
- Hydrophobic_MolePct
- Polar_minus_Nonpolar
- Amino acid composition (20 residues)
- Property group mole percentages
- Optional Dayhoff statistics (if `-d` is used)

---

## Trimmed Output Mode

The `-t` / `--trimmed` option outputs a compact feature table with:

- Gene_ID
- Residues
- Average_Residue_Weight
- Charge
- Isoelectric_Point
- Charge_per_residue
- Hydrophobic_MolePct
- Polar_minus_Nonpolar
- Tiny_MolePct
- Small_MolePct
- Aliphatic_MolePct
- Aromatic_MolePct
- Non_polar_MolePct
- Polar_MolePct
- Charged_MolePct
- Basic_MolePct
- Acidic_MolePct

---

## Quick Test

A small example file `sample.fasta` is included.

Run:

```bash
./bin/pepstats-run -f sample.fasta -o test_output/
```

Check results:

```bash
ls test_output/tables/
```

---

## Notes

- Only `pepstats` depends on `$PATH`
- All internal scripts resolve each other via relative paths
- The parser uses only Python standard library modules
- No hidden dependencies
- Output is deterministic and structured for downstream analysis

---

## Version

v1.0.0

---

## Author

Adekola Owoyemi

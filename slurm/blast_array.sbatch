#!/bin/bash
#SBATCH --job-name=phage-host-blast
#SBATCH --output=logs/blast_%A_%a.out
#SBATCH --error=logs/blast_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --array=0-19

set -euo pipefail

# Make BLAST+ and the Python environment for this repository available here.
# Cluster-specific module names/accounts/partitions are intentionally omitted.

mkdir -p logs results/blast

python scripts/run_blast.py \
  --phage-manifest work/phages.tsv \
  --db-prefix work/blast/hosts \
  --output-dir results/blast \
  --threads "${SLURM_CPUS_PER_TASK:-1}" \
  --chunk-count 20 \
  --chunk-index "${SLURM_ARRAY_TASK_ID}"

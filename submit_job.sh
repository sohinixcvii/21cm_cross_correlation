#!/bin/bash
# =============================================================================
# submit_job.sh — SLURM submission script for run_simulation.py
# =============================================================================
# Usage:  sbatch submit_job.sh
#
# Adjust the #SBATCH directives below to match your cluster's partitions,
# memory limits, and wall-time policies. The 21cmFAST lightcone at
# HII_DIM=128 typically requires ~32–64 GB RAM and several CPU hours.
# =============================================================================

#SBATCH --job-name=21cm_lightcone
#SBATCH --output=outputs/slurm_%j.out
#SBATCH --error=outputs/slurm_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=compute          # ← change to your cluster's partition name

# ── Create outputs directory if it does not exist ────────────────────────────
mkdir -p outputs

# ── Activate the 21cmfast conda environment ──────────────────────────────────
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate 21cmfast

# ── Move to the project directory ────────────────────────────────────────────
cd "$(dirname "$(realpath "$0")")" || exit 1

# ── Run the simulation ────────────────────────────────────────────────────────
echo "Starting simulation at $(date)"
python run_simulation.py
echo "Simulation finished at $(date)"

#!/bin/bash
# SLURM submission script for SPD simulation (full run)
# Submits 250 parallel jobs, each running all 27 noise combinations sequentially.
# Usage: sbatch --array=1-250 submit_spd.sh

#SBATCH --partition=main
#SBATCH --requeue
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=logs_spd_phaseonly/%x_%A_%a.out
#SBATCH --error=logs_spd_phaseonly/%x_%A_%a.err

# If using conda for R, uncomment and adjust the two lines below
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate r440
cd "$(dirname "$0")"  # cd to script directory for relative source() paths

SEED=${SLURM_ARRAY_TASK_ID}

for W in 0.0 0.5 1.0; do
  for D in 0.0 0.5 1.0; do
    for P in 0.0 0.5 1.0; do
      echo "[$(date)] seed=${SEED} W=${W} D=${D} P=${P}"
      Rscript --vanilla run_cluster_spd_phaseonly.R \
        --seed=${SEED} --sigma_warp=${W} --sigma_dist=${D} --sigma_pert=${P} \
        --result_root=results_spd_phaseonly
    done
  done
done
echo "[$(date)] DONE seed=${SEED}"

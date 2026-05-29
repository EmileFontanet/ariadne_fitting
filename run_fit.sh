#!/bin/bash
#SBATCH --job-name=ariadne
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --array=0-640%50
#SBATCH --partition=public-cpu

module load Miniconda3

# Somehow needed to activate conda env
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ariadne_env

STAR=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" stars.txt)

echo "Running fit for $STAR"
python get_star_params.py "$STAR"


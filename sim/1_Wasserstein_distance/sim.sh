#!/bin/bash 
#SBATCH --job-name=W1
#SBATCH --array 1-750
#SBATCH --output out/W1-%A_%a.out
#SBATCH --error out/W1-%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mem=2G
#SBATCH --partition=medium

module purge
module load r

R CMD BATCH --vanilla "sim.R" "out/$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID.Routput"


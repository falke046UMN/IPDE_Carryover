#!/bin/bash -l
#SBATCH --job-name=IPDEsim_prior_0.6
#SBATCH --array=1-24
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15g
#SBATCH --tmp=10g
#SBATCH --output=logs/%x_%j.out

module load R
Rscript sims.R

#!/bin/bash -l
#SBATCH --job-name=IPDEsim
#SBATCH --array=1-20
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=15g
#SBATCH --tmp=10g
cd ~/IPDEsim

module load R/4.2.2-openblas
Rscript sims.R

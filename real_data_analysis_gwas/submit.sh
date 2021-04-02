#!/bin/bash -l
#SBATCH --job-name=simu
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu

#this file is to run JWAS
module load julia
srun julia GWAS.jl
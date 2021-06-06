#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --partition=
#SBATCH --mail-type=
#SBATCH --mail-user=

#this file is to run JWAS
module load julia
srun julia GWAS.jl
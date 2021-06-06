#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --partition=
#SBATCH --mail-type=
#SBATCH --mail-user=


module load julia


srun julia /home/tianjing/paper_simu/accuracy.jl $1 $2
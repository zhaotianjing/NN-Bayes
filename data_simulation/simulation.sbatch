#!/bin/bash -l
#SBATCH --job-name=accuracy
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500
#SBATCH --time=99:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=

module load julia

srun julia /home/tianjing/paper_simu/simulation.jl $1
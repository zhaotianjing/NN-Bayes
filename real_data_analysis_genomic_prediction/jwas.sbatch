#!/bin/bash -l
#SBATCH --job-name=simu
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=99:00:00
#SBATCH --partition=bmh
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu


module load julia


srun julia /group/qtlchenggrp/tianjing/paper_real/jwas.jl $1 $2 $3


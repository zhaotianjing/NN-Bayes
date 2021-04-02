#!/bin/bash -l
#SBATCH --job-name=simu
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=99:00:00
#SBATCH --partition=bmm
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=


module load julia


srun julia /group/qtlchenggrp/tianjing/paper_real/accuracy.jl $1 $2 $3


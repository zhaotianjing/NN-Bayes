#!/bin/bash -l
#SBATCH --job-name=simu
#SBATCH --cpus-per-task=1
#SBATCH --mem=2500
#SBATCH --time=99:00:00
#SBATCH --partition=high
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=tjzhao@ucdavis.edu


module load R


srun Rscript /group/qtlchenggrp/tianjing/GWAS/roc.R $1 $2


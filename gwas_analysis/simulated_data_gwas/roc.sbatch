#!/bin/bash -l
#SBATCH --job-name=
#SBATCH --cpus-per-task=
#SBATCH --mem=
#SBATCH --time=
#SBATCH --partition=
#SBATCH --mail-type=
#SBATCH --mail-user=


module load R


srun Rscript /group/qtlchenggrp/tianjing/GWAS/roc.R $1 $2


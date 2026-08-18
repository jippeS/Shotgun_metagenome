#!/bin/bash
#SBATCH --output=taxa_barplot_%j.out
#SBATCH --job-name=taxa_barplot
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

qiime taxa barplot --i-table $1 --i-taxonomy $2 --o-visualization $3

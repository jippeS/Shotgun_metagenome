#!/bin/bash
#SBATCH --output=kraken_features%j.out
#SBATCH --job-name=kraken_features
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn1

qiime moshpit kraken2-to-features --i-reports $1 --o-table $2 --o-taxonomy $3 --verbose

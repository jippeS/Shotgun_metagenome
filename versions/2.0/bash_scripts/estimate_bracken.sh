#!/bin/bash
#SBATCH --output=estimate_bracken_%j.out
#SBATCH --job-name=estimate_bracken_
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

qiime moshpit estimate-bracken \
        --i-kraken-reports $1 \
        --i-bracken-db $2 \
        --p-read-len $3 \
        --p-level 'S' \
        --o-reports $4 \
        --o-taxonomy $5 \
        --o-table $6 \
        --verbose 

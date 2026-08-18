#!/bin/bash
#SBATCH --output=Classifying_kraken_%j.out
#SBATCH --job-name=Classifying_kraken
#SBATCH --partition=Bytesflex
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

export TMPDIR=/export/jippe/jsil/scratch/qiime_tmp
mkdir -p $TMPDIR

qiime moshpit classify-kraken2 --i-seqs $1 --i-kraken2-db $2 --p-confidence $3 --p-threads 16 --output-dir $4 --o-hits  $5 --verbose

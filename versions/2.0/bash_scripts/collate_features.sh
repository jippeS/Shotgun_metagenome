#!/bin/bash
#SBATCH --output=collate_features%j.out
#SBATCH --job-name=collate_features
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

# last argument = output
output="${@: -1}"

# all but last = inputs
inputs=("${@:1:$#-1}")

qiime moshpit collate-kraken2-reports \
    --i-kraken2-reports "${inputs[@]}" \
    --o-collated-kraken2-reports "$output"

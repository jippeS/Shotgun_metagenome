#!/bin/bash
#SBATCH --output=megahit_%j.out
#SBATCH --job-name=megahit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cn2

    
qiime assembly assemble-megahit --i-seqs $1 --p-presets "meta-sensitive" --p-num-cpu-threads $3 --o-contigs $2 --verbose

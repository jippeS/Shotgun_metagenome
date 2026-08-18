#!/bin/bash
#SBATCH --output=indexing_contigs_%j.out
#SBATCH --job-name=indexing
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cn1

qiime assembly index-contigs --i-contigs $1 --o-index $2 --p-threads $3 --verbose

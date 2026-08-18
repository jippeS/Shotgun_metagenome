#!/bin/bash
#SBATCH --output=summarize_table_%j.out
#SBATCH --job-name=summarize_table
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn1

qiime feature-table summarize --i-table $1 --o-visualization $2

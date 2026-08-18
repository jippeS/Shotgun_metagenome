#!/bin/bash
#SBATCH --output=export_tools_%j.out
#SBATCH --job-name=export_tools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

qiime tools export --input-path $1 --output-path $2

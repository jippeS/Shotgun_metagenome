#!/bin/bash
#SBATCH --output=import_hits_%j.out
#SBATCH --job-name=import_hits
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

qiime tools import --input-path $1 --output-path $2 --type "SampleData[BLAST6]"

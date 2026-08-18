#!/bin/bash
#SBATCH --output=Sorting_bam%j.out
#SBATCH --job-name=Sorting_bam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn1

samtools sort -n -o $2 $1

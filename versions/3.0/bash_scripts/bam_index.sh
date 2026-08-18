#!/bin/bash
#SBATCH --output=Sorting_bam_%j.out
#SBATCH --job-name=Sorting_bam
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

samtools index -@ $1 -o $2 $3

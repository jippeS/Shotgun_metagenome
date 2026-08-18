#!/bin/bash
#SBATCH --output=fastqc_%j.out
#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

fastqc $1 -o $2

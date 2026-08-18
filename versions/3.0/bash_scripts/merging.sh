#!/bin/bash
#SBATCH --output=merging_%j.out
#SBATCH --job-name=merging
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2


cat $1 >> $2
cat $3 >> $4

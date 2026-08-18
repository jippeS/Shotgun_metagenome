#!/bin/bash
#SBATCH --output=removeg_v2_%j.out
#SBATCH --job-name=remove_g_v2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cn1


python3 python_scripts/remove_g_updated_v2.py $1 $2 --threads=8

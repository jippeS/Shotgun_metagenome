#!/bin/bash
#SBATCH --output=make_cache%j.out
#SBATCH --job-name=make_cache
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn1

qiime tools cache-create --cache $1;
echo "make_cache finished, adding bytes..." > $2

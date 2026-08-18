#!/bin/bash
#SBATCH --output=Classif_pathogen_%j.out
#SBATCH --job-name=Classif_pathogen
#SBATCH --nodelist=cn2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16

diamond blastx --db $1 --query $2 --outfmt 6 --threads $3 --max-target-seqs 1 --quiet -e 1e-10 --more-sensitive --block-size 8 --query-cover 70 --out $4

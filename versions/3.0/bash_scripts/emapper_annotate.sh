#!/bin/bash
#SBATCH --output=emapper_annotate_%j.out
#SBATCH --job-name=emapper_annotate_moshpit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

qiime moshpit eggnog-annotate --i-eggnog-hits $1 --i-eggnog-db $2 --p-num-cpus $3 --o-ortholog-annotations $4 --verbose

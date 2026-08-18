#!/bin/bash
#SBATCH --output=trimmomatic%j.out
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --nodelist=cn1

# is already made so no need to remade for now.

qiime tools cache-create \
	  --cache $database_cache;
	  
qiime moshpit build-kraken-db \
  --p-collection standard \
  --o-kraken2-database /export/databases/cache:kraken_standard \
  --o-bracken-database /export/databases/cache:bracken_standard \
  --verbose  

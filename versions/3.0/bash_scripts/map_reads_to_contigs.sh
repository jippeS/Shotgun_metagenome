#!/bin/bash
#SBATCH --output=mapping_contigs_%j.out
#SBATCH --job-name=mapping_contigs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cn2

qiime assembly _map-reads-to-contigs --i-index $1 --i-reads $2 --o-alignment-map $3 --p-threads $4 --p-sensitivity "sensitive" --p-seed 711 --verbose

#qiime assembly _map-reads-to-contigs --i-indexed $1 --i-reads $2 $3 --p-sensitivity "sensitive" --p-threads $4 --p-seed 711 --o-alignment-map $5 --verbose
#qiime assembly map-reads-to-contigs --i-indexed-contigs "${merged_dir}/cache:contigs_indexed" --i-reads "${r1} ${r2}" --p-sensitivity "sensitive" --p-threads "$cpus" --p-seed 711 --o-alignment-map "${merged_dir}/cache:read_aln_${sample}" --verbose

#!/bin/bash
#SBATCH --output=htseq_count_%j.out
#SBATCH --job-name=htseq_count
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

htseq-count -f bam -r name -s no -t CDS -i ID $1 $2 > $3

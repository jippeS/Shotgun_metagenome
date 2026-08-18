#!/bin/bash
#SBATCH --output=make_artefact%j.out
#SBATCH --job-name=make_artefact
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path $1 --input-format PairedEndFastqManifestPhred33 --output-path $2

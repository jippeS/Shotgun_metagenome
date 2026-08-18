#!/bin/bash
#SBATCH --output=trimmomatic%j.out
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --nodelist=cn1


java -jar /export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 9 $1 $2 $3 $4 $5 $6 ILLUMINACLIP:/export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/adapters/Adapters_Macrogen.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:30 HEADCROP:1 2> $7

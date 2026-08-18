#!/bin/bash
#SBATCH --output=emapper_annotate_manual_%j.out
#SBATCH --job-name=emapper_annotate_manual
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

python3 eggnog-mapper/emapper.py -i $1 -o $2 -m diamond --dmnd_db "/export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d787/data/ref_db.dmnd" --data_dir "/export/databases/cache/data/9fc589c8-fbc9-40b6-8660-a11fcedc1f50/data/" --itype metagenome --output_dir $3 --temp_dir $4 --cpu $5 --override
       

#!/bin/bash
#SBATCH --output=emapper_%j.out
#SBATCH --job-name=emap_eggnog
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8

emapper.py -i $1 -o "merged_dna" -m diamond --no_annot --dmnd_db "/export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d787/data/ref_db.dmnd" --itype metagenome --output_dir $2 --temp_dir $3 --cpu $4 --override
#/export/jippe/jsil/programs/Shotgun/versions/3.0/eggnog-mapper/eggnogmapper/
echo "Emapper finished and also finished the file" > $5



#----------------------------------------------------------------------------
## FUNCTIONAL CLASSIFICATION 

# Straight using emapper.py in order to output .gff files

# find sample-specific contigs in cache

# locate merged contigs fasta file
fa_file=$(find "${merged_dir}/cache" -type f -path "*/data/*.fa" -print -quit)

if [ -n "$fa_file" ]; then
merged_input_contigs="$fa_file"

# create outdir
mkdir -p "${merged_dir}/results/emapper_output"

# run eggnog mapper using emapper.py

## CONSULTED emapper.py call from Moshpit_Q2_DNA_109.out which was successful: Note that emapper.py calls diamond blastx internally

# emapper.py  -i "$sample_dir/$sample/cache/377895b5-8506-42ff-bd27-c63402231687/data/IDIN19_contigs.fa" -o IDIN19 -m diamond --no_
# annot --dmnd_db /export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d787/data/ref_db.dmnd --itype metagenome --output_dir /export/microlab/temp/te
# mp_109/tmpvbcqibqp --cpu "$cpus"
# /export/microlab/miniconda3/envs/qiime2-metagenome-2024.5/bin/diamond blastx -d '/export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d7
# 87/data/ref_db.dmnd' -q '/export/microlab/temp/temp_109/qiime2/microlab/data/377895b5-8506-42ff-bd27-c63402231687/data/IDIN19_contigs.fa' --threads 16 -o '
# /export/microlab/temp/temp_109/tmpvbcqibqp/IDIN19.emapper.hits' --tmpdir '/export/microlab/users/PVEE/IDIN/DNA_trials/emappertmp_dmdn_fdcigtt3' --sensitive
#  --iterate -e 0.001 --max-target-seqs 0 --max-hsps 0  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovh
# sp scovhsp

	if [[ ! -f "${merged_dir}/results/emapper_output/merged_dna.emapper.genepred.gff" ]]; then
	# run emapper.py if output gff file does not exist  
	emapper.py  \
	  -i "$merged_input_contigs" \
	  -o "merged_dna" \
	  -m diamond \
	  --no_annot \
	  --dmnd_db "/export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d787/data/ref_db.dmnd" \
	  --itype metagenome \
	  --output_dir "${merged_dir}/results/emapper_output" \
	  --temp_dir "/export/microlab/temp/temp_${SLURM_JOB_ID}" \
	  --cpu "$cpus"

		echo -e "Eggnog output from emapper.py stored in ${merged_dir}/results/emapper_output" 
	else
		echo -e "GFF gene prediction file already exists for $sample"
	fi 

else
	echo "No .fa file found in the cache directory for merged data"
	echo "Functional classification of contigs not executed using emapper.py for merged data"
fi



#----------------------------------------------------------------------------
# MERGE READS ACROSS SAMPLES

# merge all reads.qza files from samples in sample_list into merged_reads.qza

# Define output files for the combined FASTQ files
combined_R1="${merged_dir}/combined_R1.fastq.gz"
combined_R2="${merged_dir}/combined_R2.fastq.gz"

# Concatenate the forward reads
if [[ ! -f "$combined_R1" ]]; then
	echo "combining reads from multiple samples..."
	for sample in "${sample_list[@]}"; do
		# Define the directory where the reads are stored
		reads_file1=$(find "${sample_dir}/${sample}/cache" -type f -path "*/data/*R1_001.fastq.gz" -print -quit)
	    cat "${reads_file1}" >> "$combined_R1"
	done

	# Concatenate the reverse reads
	for sample in "${sample_list[@]}"; do
		reads_file2=$(find "${sample_dir}/${sample}/cache" -type f -path "*/data/*R2_001.fastq.gz" -print -quit)
	    cat "${reads_file2}" >> "$combined_R2"
	done

	echo "Combined FASTQ files created: $combined_R1 and $combined_R2"
else
	echo "combined_reads .qza files already exist"
fi

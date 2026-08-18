
#----------------------------------------------------------------------------
## COUNT GENE ABUNDANCE

 # Sorted alignment file
sorted_alignment="${merged_dir}/results/emapper_output/$merged_dna_aln.sorted.bam"

if [[ ! -f "$sorted_alignment" ]]; then
    # Locate and copy the BAM file to emapper_results directory
    bam_file_orig=$(find "${merged_dir}/cache" -type f -path "*/data/*.bam" -print -quit)
    if [[ -n "$bam_file_orig" ]]; then
        bam_file="${merged_dir}/results/emapper_output/$merged_dna_alignment.bam"
        cp "$bam_file_orig" "$bam_file"
        echo "BAM file located and copied for merged data"

        # Sort the BAM file by read name and create the sorted alignment
        echo "sorting BAM file by read name..."
        samtools sort -n -o "$sorted_alignment" "$bam_file"
    
        # Index the copied BAM file
        echo "Indexing sorted BAM file..."
        samtools index "$sorted_alignment" -@ "$cpus"

        echo "BAM file copied, sorted by read name (-n) and indexed for sample $sample"

    else
        echo "No BAM file found for sample $sample"
    fi
else
    echo "A sorted alignment already exists for merged data: $sorted_alignment"
fi

echo "BAM file ready for counting read abundance of functional genes using HTseq-count for merged data"



# Try to deactivate the environment using the available command
if deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'deactivate'."
elif conda deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'conda deactivate'."
elif source deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'source deactivate'."
else
    echo "No deactivate command succeeded."
fi

sleep 2

# activate base with htseq
echo "activating base conda for htseq-count function..."
conda activate base

# go to merged data folder
cd "${merged_dir}/results/emapper_output"

# run htseq-count 
echo "Run HTseq-count on aln.sorted.bam and emapper.py output emapper.genepred.gff file"

# define output gene abundance count file
counts="${merged_dir}/results/emapper_output/$merged_dna_HTSeq_gene_counts.txt"

if [ ! -f "$counts" ]; then

    echo "running htseq-count for sample $sample"

    # Capture start time
    start_time=$SECONDS

    # define files
    sorted_alignment="${merged_dir}/results/emapper_output/merged_dna_aln.sorted.bam"
    gff_file="${merged_dir}/results/emapper_output/$merged_dna.emapper.genepred.gff"
    
    # run htseq-count
    htseq-count -f bam -r name -s no -t CDS -i ID "$sorted_alignment" "$gff_file" > "$counts"

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime counting read abudance of functional genes for merged data: $((runtime / 60)) minutes and $((runtime % 60)) seconds"

else
    echo "Gene abundance results already exist for merged data: $counts"
fi



# Try to deactivate the environment using the available command
if deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'deactivate'."
elif conda deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'conda deactivate'."
elif source deactivate 2>/dev/null; then
    echo "Successfully deactivated using 'source deactivate'."
else
    echo "No deactivate command succeeded."
fi

# Capture process end time
end_script=$SECONDS
# Calculate runtime
total_runtime=$((end_script - start_script))
# Output runtime
echo "Total process runtime: $((total_runtime / 60)) minutes and $((total_runtime % 60)) seconds"

sleep 2


# fin

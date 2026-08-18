
#----------------------------------------------------------------------------
# ASSEMBLE READS TO CONTIGS WITH MEGAHIT


# assembly clean reads into contigs using Megahit v2
merged_contigs_key="$sample_dir/$sample/cache/keys/contigs"

if [ ! -f "$merged_contigs_key" ]; then
    echo "Assemble reads into contigs using Megahit V2 for merged_reads.qza ..."
    # Capture start time
    start_time=$SECONDS

    qiime assembly assemble-megahit \
        --i-seqs "${merged_dir}/cache:reads" \
        --p-presets "meta-sensitive" \
        --p-num-cpu-threads "$cpus" \
        --o-contigs "${merged_dir}/cache:contigs" \
        --verbose

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime assembly: $((runtime / 60)) minutes and $((runtime % 60)) seconds"    

else
    echo "Key 'contigs' already exists in cache for merged data"
fi   



#----------------------------------------------------------------------------
# INDEX CONTIGS

# Capture start time
merged_contigs_idx_key="$sample_dir/$sample/cache/keys/contigs_indexed"

if [ ! -f "$merged_contigs_idx_key" ]; then
    echo "Indexing assembled contigs for $sample..."
    # Capture start time
    start_time=$SECONDS

    # index contigs using bowtie2
    qiime assembly index-contigs \
        --i-contigs "${merged_dir}/cache:contigs" \
        --p-threads "$cpus" \
        --o-index "${merged_dir}/cache:contigs_indexed" \
        --verbose

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime of indexing contigs: $((runtime / 60)) minutes and $((runtime % 60)) seconds"

else
    echo "Key 'contigs_indexed' already exists in cache for merged data"
fi    

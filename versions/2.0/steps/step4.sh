
#----------------------------------------------------------------------------
# STORE READS IN CACHE

# store paired-end-demux.qza in cache per sample
# DNA
merged_reads_key="${merged_dir}/cache/reads"
    
if [ ! -f "$merged_reads_key" ]; then
    echo -e "Storing merged reads into cache..."
    # Capture start time
    start_time=$SECONDS

    qiime tools cache-store \
    	--cache "${merged_dir}/cache" \
    	--artifact-path "${merged_dir}/combined_reads.qza" \
    	--key reads

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime importing reads for all samples: $((runtime / 60)) minutes and $((runtime % 60)) seconds"

else
    echo "Key 'reads' already exists in cache for merged data"
fi

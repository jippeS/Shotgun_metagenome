
#----------------------------------------------------------------------------
# MAP READS TO CONTIGS 

merged_aln_key="$sample_dir/$sample/cache/keys/read_aln_contigs"

if [ ! -f "$merged_aln_key" ]; then
    echo "Map reads to contigs with MetaBat2 for merged data ..."
    # Capture start time
    start_time=$SECONDS

    # map reads to contigs bowtie2 (trim first 5 bases from reads if not done by kneaddata)
    qiime assembly map-reads-to-contigs \
        --i-indexed-contigs "${merged_dir}/cache:contigs_indexed" \
        --i-reads "${merged_dir}/cache:reads" \
        --p-sensitivity "sensitive" \
        --p-threads "$cpus" \
        --p-seed 711 \
        --o-alignment-map "${merged_dir}/cache:read_aln_contigs" \
        --verbose

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime of mapping reads to contigs: $((runtime / 60)) minutes and $((runtime % 60)) seconds"

else
    echo "Key 'read_aln_contigs' already exists in cache for merged data"
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

sleep 2

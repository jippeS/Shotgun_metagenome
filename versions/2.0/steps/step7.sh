

#----------------------------------------------------------------------------
# EVALUATE CONTIGS WITH QUAST

# # Capture start time
# start_time=$SECONDS

# mkdir -p "${merged_dir}/results"
# evaluate assemblies (multithreading only on linux # --p-threads 6 \ )
## DOES NOT WORK COMPLETELY due to internal error
    # qiime assembly evaluate-contigs \
    #     --i-contigs "${merged_dir}/cache:contigs" \
    #     --p-memory-efficient \
    #     --p-k-mer-stats \
    #     --p-threads "$cpus" \
    #     --o-visualization "${merged_dir}/results/contigs.qzv" \
    #     --verbose
        ## DOES NOT WORK COMPLETELY

# # Capture end time
# end_time=$SECONDS
# # Calculate runtime
# runtime=$((end_time - start_time))
# # Output runtime
# echo "Runtime read mapping to contigs: $((runtime_seconds / 60)) minutes and $((runtime_seconds % 60)) seconds"


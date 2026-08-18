

#----------------------------------------------------------------------------
# CREATE MANIFEST MERGED DATA

# Create the manifest file
manifest_file="${merged_dir}/combined_manifest.tsv"

# Write the manifest content
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > $manifest_file
echo -e "combined_sample\t${combined_R1}\t${combined_R2}" >> $manifest_file

#----------------------------------------------------------------------------
# IMPORT READS USING MANIFEST MERGED DATA


output_file="${merged_dir}/combined_reads.qza"
    
    if [ ! -f "$output_file" ]; then
        echo "Running qiime tools import for $sample..."
        # Capture start time
        start_time=$SECONDS

        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "${merged_dir}/combined_manifest.tsv" \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path "$output_file"

        # Capture end time
        end_time=$SECONDS
        # Calculate runtime
        runtime=$((end_time - start_time))
        # Output runtime
        echo "Runtime importing reads for combined reads: $((runtime / 60)) minutes and $((runtime % 60)) seconds"

    else
        echo -e "combined reads.qza already exists: $output_file"
    fi


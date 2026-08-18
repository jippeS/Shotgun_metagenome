

#----------------------------------------------------------------------------
## FUNCTIONAL ANNOTATION OF EMAPPER HITS 

# import annotation hits as qza

# define database
eggnog_annot="/export/databases/cache:eggnog_annotation"

# define EggNOG hits files
emapper_hits="${merged_dir}/results/emapper_output/merged_dna.emapper.hits"
emapper_hits_qza="${merged_dir}/results/emapper_output/merged_dna.emapper.hits.qza"
annotation_file="${merged_dir}/results/emapper_output/merged_dna_annotations_contigs.qza"

if [ ! -f "$emapper_hits_qza" ]; then
	# import emapper.py hits output as qza	
	qiime tools import \
	--input-path "$emapper_hits" \
	--output-path "$emapper_hits_qza" \
	--type "SampleData[BLAST6]"

	echo "imported emapper.py hits as qza file for merged data"
else
	echo "emapper hits qza already exists: $emapper_hits_qza"
fi


if [ ! -f "$annotation_file" ]; then
    echo "Annotate merged contig hits with EggNOG Gene Ontology DB..."
    # Capture start time
    start_time=$SECONDS

    qiime moshpit eggnog-annotate \
        --i-eggnog-hits "$emapper_hits_qza" \
        --i-eggnog-db "$eggnog_annot" \
        --p-num-cpus "$cpus" \
        --o-ortholog-annotations "$annotation_file" \
        --verbose

    # Capture end time
    end_time=$SECONDS
    # Calculate runtime
    runtime=$((end_time - start_time))
    # Output runtime
    echo "Runtime EggNOG annotation of merged DNA contig hits: $((runtime / 60)) minutes and $((runtime % 60)) seconds"
else
    echo "EggNOG annotations for merged DNA contig hits already exists: $annotation_file"
fi   

    # Export EggNOG annotations from the qza file
    qiime tools export \
        --input-path "$annotation_file" \
        --output-path "${merged_dir}/results/emapper_output"

    echo "Exported EggNOG annotations to ${merged_dir}/results/emapper_output"

sleep 2

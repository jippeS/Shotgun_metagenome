#!/bin/bash

inputdir=$1
index_dir=$2
index_manifest=$3

sample_list=()

for d in "$inputdir"/*/; do
    d=${d%/}
    name=$(basename "$d")

    # skip non-sample directories
    if [[ "$name" == *"_merged" ]] || \
       [[ "$name" == temp* ]] || \
       [[ "$name" == cache ]] || \
       [[ "$name" == output ]] || \
       [[ "$name" == results ]]; then
        continue
    fi

    sample_list+=("$name")
done

echo -e "sample-id\tindex-path" > "$index_manifest"

for sample in "${sample_list[@]}"; do
    echo -e "${sample}\t${index_dir}" >> "$index_manifest"
done

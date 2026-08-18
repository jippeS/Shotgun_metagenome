#!/bin/bash

inputdir=$1
fw_reads=$2
rv_reads=$3
manifest=$4

# correct QIIME2 header
echo -e "sample-id,absolute-filepath,direction" > "$manifest"

for d in "$inputdir"/*/; do
    d=${d%/}
    sample=$(basename "$d")

    # skip non-sample dirs
    if [[ "$sample" == *"_merged" ]] || \
       [[ "$sample" == temp* ]] || \
       [[ "$sample" == cache ]] || \
       [[ "$sample" == output ]] || \
       [[ "$sample" == results ]]; then
        continue
    fi

    # forward
    echo -e "${sample},${fw_reads},forward" >> "$manifest"

    # reverse
    echo -e "${sample},${rv_reads},reverse" >> "$manifest"
done

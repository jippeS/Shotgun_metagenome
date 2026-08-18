configfile: "config.yaml"

from os import listdir
from os.path import isfile, join
import os
import re
import argparse
import sys

# Make the outputdir
startdir= config["startdir"]

def create_dirs(key_dir, key_dir_output):
    if not os.path.exists(key_dir):
        os.makedirs(key_dir)
    if not os.path.exists(key_dir_output):
        os.makedirs(key_dir_output)

def organize_files_by_prefix(startdir):

    # Get the list of file names in the specified folder
    file_names = [name for name in os.listdir(startdir) if os.path.isfile(os.path.join(startdir, name)) and name.endswith(".fastq.gz")]
    input_wildcard = [name.split("_")[0] for name in file_names]
    input_wildcard = list(set(input_wildcard))
    # Populate dic_dir_with_files with files from the inputdir
    for file in file_names:
        key = file.split("_")[0]
        # Make the directories
        key_dir = os.path.join(startdir, key, "input/")
        key_dir_output = os.path.join(startdir, key, "output/")

        # create benchmark folder
        path = key_dir_output + "benchmarks/Importing_data.txt"
        os.makedirs(path,exist_ok=True)

        #create input and output dirs
        create_dirs(key_dir,key_dir_output)

        forward_or_reverse = file.split("_")[-1]
        if forward_or_reverse.startswith("1"):
            file1 = "forward.fastq.gz"

        elif forward_or_reverse.startswith("2"):
            file1 = "reverse.fastq.gz"
        else:
            sys.exit(print("Look at code def organize_file_by_prefix to fix renaming the headers."))
        new_file_destination = os.path.join(key_dir, os.path.basename(file1))
        old_path_file = os.path.join(startdir, file)
        # Move the file to the new folder
        try:
            os.rename(old_path_file, new_file_destination)
        except:
            continue

    # check for already made directories which are not called output.
    if len(input_wildcard) < 1:
        input_wildcard = [name for name in os.listdir(startdir) if os.path.isdir(os.path.join(startdir, name)) and name != "output" and name != ".snakemake" and name != "multiqc" and name != "post_multiqc"]

    return input_wildcard



input_wildcard = organize_files_by_prefix(startdir)
print(input_wildcard)
# input_wildcard = ['mGMAC14']
# print(input_wildcard)

rule all:
    input:
        expand(startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.fsa", dataset=input_wildcard),
        expand(startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.fsa", dataset=input_wildcard),
        expand(startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.fsa", dataset=input_wildcard),
        expand(startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.fsa", dataset=input_wildcard),
        expand(startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.fsa", dataset=input_wildcard)

rule fastqc_forward:
    output:
        forward_fastqc = startdir + "{dataset}/output/fastqc/forward_fastqc.html",
        forward_fastqc_zip = startdir + "{dataset}/output/fastqc/forward_fastqc.zip"
    conda:
        "multiqc.yaml"
    params:
        forward_file =  startdir + "{dataset}/input/forward.fastq.gz",
        output_dir = startdir + "{dataset}/output/fastqc/"
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/forward.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.forward_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'

rule fastqc_reverse:
    input:
        fw_fastqc = rules.fastqc_forward.output.forward_fastqc,
        fw_fastqc_zip = rules.fastqc_forward.output.forward_fastqc_zip
    output:
        reverse_fastqc = startdir + "{dataset}/output/fastqc/reverse_fastqc.html",
        reverse_fastqc_zip = startdir + "{dataset}/output/fastqc/reverse_fastqc.zip"
    conda:
        "multiqc.yaml"
    params:
        reverse_file = startdir + "{dataset}/input/reverse.fastq.gz",
        output_dir = startdir + "{dataset}/output/fastqc/",
        done_folder = startdir +"{dataset}/output/done/"
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/reverse.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.reverse_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'
       'mkdir {params.done_folder}'

rule multiqc:
    input:
        expand(startdir + "{dataset}/output/fastqc/reverse_fastqc.zip", dataset=input_wildcard)
    output:
        startdir + "multiqc/multiqc_report.html"
    conda:
        "multiqc.yaml"
    params:
        outputdir = startdir + "multiqc/"
    message:
        "@#"
        "Multiqc:   "
        "multiqc "
        "   */output/fastqc/* "
        "   --outdir multiqc/ "
        "   -f "
        "   -d "
        "   -s"
        "@#"
    shell:
        "multiqc {startdir}*/output/fastqc/* --outdir {params.outputdir} -f -d -s;"

rule pre_process:
    input:
        rules.multiqc.output
    output:
        touch(startdir + "{dataset}/output/done/1_pre_process_done.txt")
    conda:
        "multiqc.yaml"
    params:
        output_folder = startdir + "{dataset}/output/processed/",
        file_path1 = startdir + "{dataset}/input/forward.fastq.gz",
        forward_file = startdir + "{dataset}/output/processed/forward_filtered.fastq.gz"
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/forward.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path1} {params.output_folder};"
        "python3 python_scripts/remove_g_updated.py {params.file_path1} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {params.forward_file} --seconds=300;"

rule pre_process2:
    input:
        rules.pre_process.output
    output:
        touch(startdir + "{dataset}/output/done/2_pre_process_done.txt")
    conda:
        "multiqc.yaml"
    params:
        output_folder = startdir + "{dataset}/output/processed/",
        file_path2 = startdir + "{dataset}/input/reverse.fastq.gz",
        reverse_file = startdir + "{dataset}/output/processed/reverse_filtered.fastq.gz"
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/reverse.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path2} {params.output_folder};"
        "python3 python_scripts/remove_g_updated.py {params.file_path2} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {params.reverse_file} --seconds=300;"

rule trimmomatic:
    input:
        rules.pre_process2.output
    output:
        output_log = startdir + "{dataset}/output/processed/trimmomatic.log",
        output_2 = startdir + "{dataset}/output/done/3_trimmomatic_done.txt"
    conda:
        "multiqc.yaml"
    params:
        forward_file=startdir + "{dataset}/output/processed/forward_filtered.fastq.gz",
        reverse_file=startdir + "{dataset}/output/processed/reverse_filtered.fastq.gz",
        forward_output = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz",
        reverse_output = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz",
        forward_output_U = startdir + "{dataset}/output/processed/forward_filtered_U.fastq.gz",
        reverse_output_U = startdir + "{dataset}/output/processed/reverse_filtered_U.fastq.gz"
    message:
        "@#"
        "Trimming:   "
        "java -jar /export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar   "
        "   PE"
        "   -threads 1"
        "   */output/processed/forward_filtered.fastq.gz */output/processed/reverse_filtered.fastq.gz */output/processed/forward_filtered_P.fastq.gz */output/processed/forward_filtered_U.fastq.gz */output/processed/reverse_filtered_P.fastq.gz */output/processed/reverse_filtered_U.fastq.gz"
        "   ILLUMINACLIP:/export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/adapters/Adapters_Macrogen.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:30 HEADCROP:1 2> */output/processed/trimmomatic.log"
        "@#"
    shell:
       'sbatch bash_scripts/trimmomatic.sh {params.forward_file} {params.reverse_file} {startdir}{wildcards.dataset}/output/processed {output.output_log};'
       'python3 {config[tooldir]}wait_file.py {params.forward_output} {params.reverse_output} --seconds=200;'
       'rm {params.forward_output_U}; rm {params.reverse_output_U};'
       'rm {params.forward_file}; rm {params.reverse_file};'
       'touch {output.output_2}'

rule post_fastqc_forward:
    input:
        trimmomatic = rules.trimmomatic.output.output_log
    output:
        forward_fastqc = startdir + "{dataset}/output/post_fastqc/forward_filtered_P_fastqc.html",
        forward_fastqc_zip = startdir + "{dataset}/output/post_fastqc/forward_filtered_P_fastqc.zip"
    conda:
        "multiqc.yaml"
    params:
        output_dir = startdir + "{dataset}/output/post_fastqc/",
        forward_file = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz"
    message:
        "@#"
        "Post trimmomatic Quality control forward:   "
        "fastqc "
        "   */output/processed/forward_filtered_P.fastq.gz "
        "   -o */output/post_fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.forward_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'

rule post_fastqc_reverse:
    input:
        fw_fastqc = rules.post_fastqc_forward.output.forward_fastqc,
        fw_fastqc_zip = rules.post_fastqc_forward.output.forward_fastqc_zip
    output:
        reverse_fastqc = startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.html",
        reverse_fastqc_zip = startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.zip"
    conda:
        "multiqc.yaml"
    params:
        output_dir = startdir + "{dataset}/output/post_fastqc/",
        reverse_file = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz"
    message:
        "@#"
        "Post trimmomatic Quality control reverse:   "
        "fastqc "
        "   */output/processed/reverse_filtered_P.fastq.gz "
        "   -o */output/post_fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.reverse_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'

rule post_multiqc:
    input:
        elly = expand(startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.zip", dataset=input_wildcard)
    output:
        startdir + "post_multiqc/multiqc_report.html"
    conda:
        "multiqc.yaml"
    params:
        outputdir = startdir + "post_multiqc/"
    message:
        "@#"
        "Post processing Multiqc:   "
        "multiqc "
        "   */output/fastqc/* "
        "   --outdir post_multiqc/ "
        "   -f "
        "   -d "
        "   -s"
        "@#"
    shell:
         "multiqc {startdir}*/output/post_fastqc/* --outdir {params.outputdir} -f -d -s"

rule bbmap_interleave:
    input:
        input1 = rules.trimmomatic.output.output_log,
        input2 = rules.post_multiqc.output
    output:
        startdir + "{dataset}/output/processed/interleaved.fastq.gz"
    conda:
        "multiqc.yaml"
    params:
        forward_output = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz",
        reverse_output = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz",
    message:
        "@#"
        "Interleaving:   "
        "reformat.sh"
        "   in1=*/output/processed/forward_filtered_P.fastq.gz"
        "   in2=*/output/processed/reverse_filtered_P.fastq.gz"
        "   out=*/output/processed/interleaved.fastq.gz"
        "@#"
    shell:
       'sbatch bash_scripts/bbmap_interleave.sh {params.forward_output} {params.reverse_output} {output};'
       'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'
       'rm {params.forward_output}; rm {params.reverse_output}'

rule classify_kma:
    input:
        rules.bbmap_interleave.output
    output:
        output1 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.fsa",
        output2 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.frag.gz",
        output3 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.aln",
        output4 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.res",
        output5 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.tsv"
    conda:
        "multiqc.yaml"
    params:
        output_name = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output"
    message:
        "@#"
        "Classifying card:   "
        "kma "
        "   -int $1"
        "   -o $2"
        "   -t_db /export/databases/Wetsus/Card/3.2.7/kma_index/nucleotide_fasta_protein_homolog_model_variants"
        "   -tsv"
        "   -1t1"
        "   -cge"
        "   -hmm"
        "   -ex_mode"
        "   -mem_mode"
        "@#"
    shell:
        'sbatch bash_scripts/kma_classify.sh {input} {params.output_name};'
        'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'

rule classify_kma_resfinder:
    input:
        rules.bbmap_interleave.output
    output:
        output1 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.fsa",
        output2 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.frag.gz",
        output3 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.aln",
        output4 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.res",
        output5 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.tsv"
    conda:
        "multiqc.yaml"
    params:
        output_name = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output"
    message:
        "@#"
        "Classifying Resfinder:   "
        "kma "
        "   -int $1"
        "   -o $2"
        "   -t_db /export/databases/Wetsus/Resfinder2/ResFinder2"
        "   -tsv"
        "   -1t1"
        "   -cge"
        "   -hmm"
        "   -ex_mode"
        "   -mem_mode"
        "@#"
    shell:
        'sbatch bash_scripts/kma_classify_resfinder.sh {input} {params.output_name};'
        'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'

rule classify_kma_silva:
    input:
        rules.bbmap_interleave.output
    output:
        output1 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.fsa",
        output2 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.frag.gz",
        output3 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.aln",
        output4 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.res",
        output5 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.tsv"
    conda:
        "multiqc.yaml"
    params:
        output_name = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output"
    message:
        "@#"
        "Classifying Resfinder:   "
        "kma "
        "   -int $1"
        "   -o $2"
        "   -t_db /export/databases/Silva/138.2/kma/silva_138.2_nr99_kma"
        "   -tsv"
        "   -1t1"
        "   -cge"
        "   -hmm"
        "   -ex_mode"
        "   -mem_mode"
        "@#"
    shell:
        'sbatch bash_scripts/kma_classify_silva138.2.sh {input} {params.output_name};'
        'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'

rule classify_kma_panres:
    input:
        rules.bbmap_interleave.output
    output:
        output1=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.fsa",
        output2=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.frag.gz",
        output3=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.aln",
        output4=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.res",
        output5=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.tsv"
    conda:
        "multiqc.yaml"
    params:
        output_name=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output"
    message:
        "@#"
    "Classifying Resfinder:   "
    "kma "
    "   -int $1"
    "   -o $2"
    "   -t_db /export/databases/Wetsus/Panres/kma/panres_31_07_2025"
    "   -tsv"
    "   -1t1"
    "   -cge"
    "   -hmm"
    "   -ex_mode"
    "   -mem_mode"
    "@#"
    shell:
        'sbatch bash_scripts/kma_panres.sh {input} {params.output_name};'
        'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'


rule classify_kma_pr2:
    input:
        rules.bbmap_interleave.output
    output:
        output1=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.fsa",
        output2=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.frag.gz",
        output3=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.aln",
        output4=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.res",
        output5=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.tsv"
    conda:
        "multiqc.yaml"
    params:
        output_name=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output"
    message:
        "@#"
    "Classifying Resfinder:   "
    "kma "
    "   -int $1"
    "   -o $2"
    "   -t_db /export/databases/Wetsus/PR2/5.1.0/kma_16/pr2_5.1.0_SSU"
    "   -tsv"
    "   -1t1"
    "   -cge"
    "   -hmm"
    "   -ex_mode"
    "   -mem_mode"
    "@#"
    shell:
        'sbatch bash_scripts/kma_pr2.sh {input} {params.output_name};'
        'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
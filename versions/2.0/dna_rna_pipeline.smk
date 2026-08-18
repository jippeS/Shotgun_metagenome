configfile: "config.yaml"

from os import listdir
from os.path import isfile, join
import os
import re
import argparse
import sys

folder_path_DNA = config["inputdir"] + "DNA/"
folder_path_RNA = config["inputdir"] + "DNA/"

ignore = {"output", "logs", ".snakemake", "temp"}

#Retrieve the files.
dna_files = [f for f in os.listdir(folder_path_DNA) if os.path.isdir(os.path.join(folder_path_DNA, f)) and f not in ignore]
rna_files = [f for f in os.listdir(folder_path_DNA) if os.path.isdir(os.path.join(folder_path_DNA, f)) and f not in ignore]


conf_thresh = str(config["CONF_THRESH"] )

rule all:
    input:
        expand(config["outputdir"] + "DNA/{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza", fq_file=dna_files),
        expand(config["outputdir"] + "RNA/{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza", fq_file=rna_files),
        # expand(config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_table.qzv", conf_thresh=conf_thresh),
        # expand(config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_barplot.qzv", conf_thresh=conf_thresh)

rule fastqc_forward_dna:
    output:
        forward_fastqc = config["outputdir"] + "DNA/pre_fastqc/{fq_file}_1_fastqc.html",
        forward_fastqc_zip = config["outputdir"] + "DNA/pre_fastqc/{fq_file}_1_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir = config["outputdir"] + "DNA/pre_fastqc/",
        input_file = config["inputdir"] + "DNA/{fq_file}/input/forward.fastq.gz"
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/forward.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.input_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'
#
rule fastqc_reverse_dna:
    input:
        input_file2 = rules.fastqc_forward_dna.output.forward_fastqc_zip
    output:
        reverse_fastqc = config["outputdir"] + "DNA/output/pre_fastqc/{fq_file}_2_fastqc.html",
        reverse_fastqc_zip = config["outputdir"] + "DNA/output/pre_fastqc/{fq_file}_2_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir = config["outputdir"] + "DNA/output/pre_fastqc/",
        input_file = config["inputdir"] + "DNA/{fq_file}/input/reverse.fastq.gz"
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/reverse.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
       'sbatch bash_scripts/fastq.sh {params.input_file} {params.output_dir};'
       'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'

rule multiqc_dna:
    input:
        expand(config["outputdir"] + "DNA/output/pre_fastqc/{fq_file}_2_fastqc.zip", fq_file=fq_files)
    output:
        config["outputdir"] + "DNA/output/pre_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir = config["outputdir"] + "DNA/output/pre_multiqc/",
        inputtdir = config["outputdir"] + "DNA/output/pre_fastqc/"
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
        "multiqc {params.inputtdir} --outdir {params.outputdir} -f -d -s"


rule pre_process1_dna:
    input:
        rules.multiqc_dna.output
    output:
        config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_1_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder = config["outputdir"] + "DNA/{fq_file}/output/processed/",
        file_path1 = config["outputdir"] + "DNA/{fq_file}/input/{fq_file}_1.fastq.gz",
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/forward.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path1} {params.output_folder};"
        "sbatch bash_scripts/remove_g_v2.sh {params.file_path1} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=300;"

rule pre_process2_dna:
    input:
        rules.pre_process1_dna.output
    output:
        config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_2_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder = config["outputdir"] + "DNA/{fq_file}/output/processed/",
        file_path1 = config["outputdir"] + "DNA/{fq_file}/input/{fq_file}_2.fastq.gz",
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/reverse.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path1} {params.output_folder};"
        "sbatch bash_scripts/remove_g_v2.sh {params.file_path1} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=300;"



rule trimmomatic_dna:
    input:
        reverse_file = rules.pre_process2.output,
        forward_file = rules.pre_process1.output
    output:
        output_log = config["outputdir"] + "DNA/{fq_file}/output/processed/trimmomatic.log",
        output_manifest = config["outputdir"] + "DNA/{fq_file}/input/pe-33v1-manifest.tsv"
    conda:
        "environments/multiqc.yaml"
    params:
        forward_output_P = config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_1_filtered_P.fastq.gz",
        reverse_output_P = config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_2_filtered_P.fastq.gz",
        forward_output_U = config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_1_filtered_U.fastq.gz",
        reverse_output_U = config["outputdir"] + "DNA/{fq_file}/output/processed/{fq_file}_2_filtered_U.fastq.gz",
        id = "{fq_file}",
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
       'sbatch bash_scripts/trimmomatic.sh {input.forward_file} {input.reverse_file} {params.forward_output_P} {params.forward_output_U}  {params.reverse_output_P} {params.reverse_output_U} {output.output_log};'
       'python3 {config[tooldir]}wait_file.py {params.forward_output_P} {params.reverse_output_P} --seconds=200;'
       'printf "sample-id,absolute-filepath,direction\n{params.id},{params.forward_output_P},forward\n{params.id},{params.reverse_output_P},reverse\n" > {output.output_manifest}'


rule post_fastqc_forward_dna:
    input:
        rules.trimmomatic_dna.output.output_manifest
    output:
        forward_fastqc = config["outputdir"] + "DNA/output/post_fastqc/{fq_file}_1_filtered_P_fastqc.html",
        forward_fastqc_zip = config["outputdir"] + "DNA/output/post_fastqc/{fq_file}_1_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "DNA/output/post_fastqc/",
        input_file=config["outputdir"] + "DNA/{fq_file}/input/{fq_file}_1.fastq.gz",
        input_file2=rules.trimmomatic_dna.params.forward_output_P
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/forward.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file2} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'

#
rule post_fastqc_reverse_dna:
    input:
        input_file2=rules.post_fastqc_forward_dna.output.forward_fastqc_zip
    output:
        reverse_fastqc=config["outputdir"] + "DNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.html",
        reverse_fastqc_zip=config["outputdir"] + "DNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "DNA/output/post_fastqc/",
        input_file=config["outputdir"] + "DNA/{fq_file}/input/{fq_file}_2.fastq.gz",
        input_file2=rules.trimmomatic_dna.params.reverse_output_P
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/reverse.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file2} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'

rule post_multiqc_dna:
    input:
        expand(config["outputdir"] + "DNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip",fq_file=fq_files)
    output:
        config["outputdir"] + "DNA/output/post_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir=config["outputdir"] + "DNA/output/post_multiqc/",
        inputtdir=config["outputdir"] + "DNA/output/post_fastqc/"
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
        "multiqc {params.inputtdir} --outdir {params.outputdir} -f -d -s"


rule make_cache_dna:
    input:
        rules.post_multiqc_dna.output
    output:
        cache_dir = directory(config["outputdir"] + "DNA/temp/{fq_file}"),
        done = config["outputdir"] + "DNA/{fq_file}/output/steps/{fq_file}_cache_made.txt"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        """
        mkdir -p {output.cache_dir}
        sbatch bash_scripts/make_cache.sh {output.cache_dir}/cache {output.done}
        python3 {config[tooldir]}wait_file.py {output.done} --seconds=40
        """

rule make_artefact_dna:
    input:
        cache = rules.make_cache_dna.output,
        manifest = rules.trimmomatic_dna.output.output_manifest
    output:
        config["outputdir"] + "DNA/{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/create_artifact.sh {input.manifest} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule fastqc_forward_rna:
    output:
        forward_fastqc = config["outputdir"] + "RNA/output/pre_fastqc/{fq_file}_1_fastqc.html",
        forward_fastqc_zip = config["outputdir"] + "RNA/output/pre_fastqc/{fq_file}_1_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "RNA/output/pre_fastqc/",
        input_file=config["outputdir"] + "{fq_file}/input/{fq_file}_1.fastq.gz"
    message:
        "@#"
    "Quality control forward:   "
    "fastqc "
    "   */input/forward.fastq.gz "
    "   -o */output/fastqc/ "
    "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file} {params.output_dir};'
    'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'

        #
rule fastqc_reverse_rna:
    input:
        input_file2=rules.fastqc_forward_rna.output.forward_fastqc_zip
    output:
        reverse_fastqc=config["outputdir"] + "RNA/output/pre_fastqc/{fq_file}_2_fastqc.html",
        reverse_fastqc_zip=config["outputdir"] + "RNA/output/pre_fastqc/{fq_file}_2_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "RNA/output/pre_fastqc/",
        input_file=config["outputdir"] + "{fq_file}/input/{fq_file}_2.fastq.gz"
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/reverse.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'

rule multiqc_rna:
    input:
        expand(config["outputdir"] + "RNA/output/pre_fastqc/{fq_file}_2_fastqc.zip",fq_file=fq_files)
    output:
        config["outputdir"] + "RNA/output/pre_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir=config["outputdir"] + "RNA/output/pre_multiqc/",
        inputtdir=config["outputdir"] + "RNA/output/pre_fastqc/"
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
        "multiqc {params.inputtdir} --outdir {params.outputdir} -f -d -s"


rule pre_process1_rna:
    input:
        rules.multiqc_rna.output
    output:
        config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_1_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder=config["outputdir"] + "RNA/{fq_file}/output/processed/",
        file_path1=config["outputdir"] + "RNA/{fq_file}/input/{fq_file}_1.fastq.gz",
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/forward.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path1} {params.output_folder};"
        "sbatch bash_scripts/remove_g_v2.sh {params.file_path1} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=300;"

rule pre_process2_rna:
    input:
        rules.pre_process1_rna.output
    output:
        config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_2_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder=config["outputdir"] + "RNA/{fq_file}/output/processed/",
        file_path1=config["outputdir"] + "RNA/{fq_file}/input/{fq_file}_2.fastq.gz",
    message:
        "@#"
        "Remove forward g:   "
        "python3 /export/jippe/jsil/programs/Transcriptoom/versions/2.0/python_scripts/remove_g_updated.py  "
        "   */input/reverse.fastq.gz"
        "   */output/processed/ &"
        "@#"
    shell:
        # "sbatch bash_scripts/remove_g_seq.sh {params.file_path1} {params.output_folder};"
        "sbatch bash_scripts/remove_g_v2.sh {params.file_path1} {params.output_folder};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=300;"


rule trimmomatic_rna:
    input:
        reverse_file=rules.pre_process2_rna.output,
        forward_file=rules.pre_process1_rna.output
    output:
        output_log=config["outputdir"] + "RNA/{fq_file}/output/processed/trimmomatic.log",
        output_manifest=config["outputdir"] + "RNA/{fq_file}/input/pe-33v1-manifest.tsv"
    conda:
        "environments/multiqc.yaml"
    params:
        forward_output_P=config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_1_filtered_P.fastq.gz",
        reverse_output_P=config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_2_filtered_P.fastq.gz",
        forward_output_U=config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_1_filtered_U.fastq.gz",
        reverse_output_U=config["outputdir"] + "RNA/{fq_file}/output/processed/{fq_file}_2_filtered_U.fastq.gz",
        id="{fq_file}",
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
        'sbatch bash_scripts/trimmomatic.sh {input.forward_file} {input.reverse_file} {params.forward_output_P} {params.forward_output_U}  {params.reverse_output_P} {params.reverse_output_U} {output.output_log};'
        'python3 {config[tooldir]}wait_file.py {params.forward_output_P} {params.reverse_output_P} --seconds=200;'
        'printf "sample-id,absolute-filepath,direction\n{params.id},{params.forward_output_P},forward\n{params.id},{params.reverse_output_P},reverse\n" > {output.output_manifest}'


rule post_fastqc_forward_rna:
    input:
        rules.trimmomatic_rna.output.output_manifest
    output:
        forward_fastqc=config["outputdir"] + "RNA/output/post_fastqc/{fq_file}_1_filtered_P_fastqc.html",
        forward_fastqc_zip=config["outputdir"] + "RNA/output/post_fastqc/{fq_file}_1_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "RNA/output/post_fastqc/",
        input_file=config["outputdir"] + "RNA/{fq_file}/input/{fq_file}_1.fastq.gz",
        input_file2=rules.trimmomatic_rna.params.forward_output_P
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/forward.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file2} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'

#
rule post_fastqc_reverse_rna:
    input:
        input_file2=rules.post_fastqc_forward_rna.output.forward_fastqc_zip
    output:
        reverse_fastqc=config["outputdir"] + "RNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.html",
        reverse_fastqc_zip=config["outputdir"] + "RNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["outputdir"] + "RNA/output/post_fastqc/",
        input_file=config["outputdir"] + "RNA/{fq_file}/input/{fq_file}_2.fastq.gz",
        input_file2=rules.trimmomatic_rna.params.reverse_output_P
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "   */input/reverse.fastq.gz "
        "   -o */output/fastqc/ "
        "@#"
    shell:
        'sbatch bash_scripts/fastq.sh {params.input_file2} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'

rule post_multiqc_rna:
    input:
        expand(config["outputdir"] + "RNA/output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip",fq_file=fq_files)
    output:
        config["outputdir"] + "RNA/output/post_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir=config["outputdir"] + "RNA/output/post_multiqc/",
        inputtdir=config["outputdir"] + "RNA/output/post_fastqc/"
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
        "multiqc {params.inputtdir} --outdir {params.outputdir} -f -d -s"


rule make_cache_rna:
    input:
        rules.post_multiqc_rna.output
    output:
        cache_dir=directory(config["outputdir"] + "RNA/temp/{fq_file}"),
        done=config["outputdir"] + "RNA/{fq_file}/output/steps/{fq_file}_cache_made.txt"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        """
            mkdir -p {output.cache_dir}
            sbatch bash_scripts/make_cache.sh {output.cache_dir}/cache {output.done}
            python3 {config[tooldir]}wait_file.py {output.done} --seconds=40
            """

rule make_artefact_rna:
    input:
        cache=rules.make_cache_rna.output,
        manifest=rules.trimmomatic_rna.output.output_manifest
    output:
        config["outputdir"] + "RNA/{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/create_artifact.sh {input.manifest} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'



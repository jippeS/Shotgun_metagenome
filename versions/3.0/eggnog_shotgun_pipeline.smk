configfile: "config.yaml"

from os import listdir
from os.path import isfile, join
import os
import re
import argparse
import sys
def is_already_organized(folder):
    """Check if folder already contains sample/input/*.fastq.gz"""
    for item in folder.iterdir():
        if item.is_dir() and item.name not in {"output", "logs", ".snakemake", "temp"}:
            input_dir = item / "input"
            if input_dir.exists():
                fastqs = list(input_dir.glob("*.fastq.gz"))
                if fastqs:
                    return True
    return False


def change_names(folder):
    mapping_file = folder / "sampleName_clientId.txt"

    suffixes_old = ["_1.fq.gz", "_2.fq.gz"]
    suffixes_new = ["_1.fastq.gz", "_2.fastq.gz"]

    dry_run = False

    if not mapping_file.exists():
        raise FileNotFoundError(f"Mapping file not found: {mapping_file}")

    with open(mapping_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            old_name, new_name = line.split("\t")

            # 📁 Create sample/input folder
            input_dir = folder / new_name / "input"
            if not dry_run:
                input_dir.mkdir(parents=True, exist_ok=True)

            for i in [0, 1]:
                old_file = folder / f"{old_name}{suffixes_old[i]}"
                new_file = input_dir / f"{new_name}{suffixes_new[i]}"

                if old_file.exists():
                    if new_file.exists():
                        print(f"SKIP (exists): {new_file}")
                        continue

                    if dry_run:
                        print(f"DRY-RUN: {old_file} → {new_file}")
                    else:
                        old_file.rename(new_file)
                        print(f"RENAMED: {old_file.name} → {new_file}")
                else:
                    print(f"MISSING: {old_file}")


def list_sample_folders(folder):
    ignore = {"output", "logs", ".snakemake", "temp"}

    return [
        f.name for f in folder.iterdir()
        if f.is_dir() and f.name not in ignore
    ]


folder = Path(config["inputdir"])

# 📁 Create output dir
outputdir = folder / "output"
outputdir.mkdir(exist_ok=True)

print(f"Input folder: {folder}")
print(f"Output folder: {outputdir}")

# 🧠 Skip if already structured
if is_already_organized(folder):
    print("\nDetected sample/input/*.fastq.gz → skipping renaming step")
else:
    print("\nOrganizing FASTQ files into sample/input/")
    change_names(folder)

# 📂 List samples
fq_files = list_sample_folders(folder)
fq_file = fq_files[0]
print("\nDetected sample folders:")
for s in fq_files:
    print(f" - {s}")


if config["CONF_THRESH"] == "DNA" or "dna":
    config["read_length"] = 150
elif config["CONF_THRESH"] == "RNA" or "rna":
    config["read_length"] = 100

conf_thresh = str(config["CONF_THRESH"] )

rule all:
    input:
        expand(config["inputdir"] + "{fq_file}/output/processed/{fq_file}_merged_HTSeq_gene_counts.txt", fq_file=fq_files)
        # expand(config["inputdir"] + "output/processed/emapper_output/merged_dna.emapper.annotations", fq_file=fq_files),

        #temporary removed because of error of indexing a namesorted bam.
        # expand(config["inputdir"] + "{fq_file}/output/processed/alignment/{fq_file}_sorted_alignment.bam.bai", fq_file=fq_files),

        # expand(config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza", fq_file=fq_files),
        # expand(config["inputdir"] + "{fq_file}/output/processed/{fq_file}_" + config["naming_convention"] + "read_aln_contigs.qza", fq_file=fq_files),

        # expand(config["inputdir"] + "{fq_file}/output/processed/emapper_output/{fq_file}_merged_dna.emapper.genepred.fasta", fq_file=fq_files),
        # expand(config["inputdir"] + "{fq_file}/output/processed/emapper_output/{fq_file}_merged_dna.emapper.genepred.gff" , fq_file=fq_files),
        # expand(config["inputdir"] + "{fq_file}/output/processed/emapper_output/{fq_file}_merged_dna_annotations_contig.qza", fq_file=fq_files),
        # config["inputdir"] + "output/processed/emapper_output/merged_dna.emapper.seed_orthologs",

        # config["inputdir"] + "output/processed/emapper_output/merged_dna_annotations_contig.qza"


rule fastqc_forward:
    output:
        forward_fastqc = config["inputdir"] + "output/pre_fastqc/{fq_file}_1_fastqc.html",
        forward_fastqc_zip = config["inputdir"] + "output/pre_fastqc/{fq_file}_1_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir = config["inputdir"] + "output/pre_fastqc/",
        input_file = config["inputdir"] + "{fq_file}/input/{fq_file}_1.fastq.gz"
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
rule fastqc_reverse:
    input:
        input_file2 = rules.fastqc_forward.output.forward_fastqc_zip
    output:
        reverse_fastqc = config["inputdir"] + "output/pre_fastqc/{fq_file}_2_fastqc.html",
        reverse_fastqc_zip = config["inputdir"] + "output/pre_fastqc/{fq_file}_2_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir = config["inputdir"] + "output/pre_fastqc/",
        input_file = config["inputdir"] + "{fq_file}/input/{fq_file}_2.fastq.gz"
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

rule multiqc:
    input:
        expand(config["inputdir"] + "output/pre_fastqc/{fq_file}_2_fastqc.zip", fq_file=fq_files)
    output:
        config["inputdir"] + "output/pre_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir = config["inputdir"] + "output/pre_multiqc/",
        inputtdir = config["inputdir"] + "output/pre_fastqc/"
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


rule pre_process1:
    input:
        rules.multiqc.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/{fq_file}_1_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder = config["inputdir"] + "{fq_file}/output/processed/",
        file_path1 = config["inputdir"] + "{fq_file}/input/{fq_file}_1.fastq.gz",
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

rule pre_process2:
    input:
        rules.pre_process1.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/{fq_file}_2_filtered.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        output_folder = config["inputdir"] + "{fq_file}/output/processed/",
        file_path1 = config["inputdir"] + "{fq_file}/input/{fq_file}_2.fastq.gz",
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



rule trimmomatic:
    input:
        reverse_file = rules.pre_process2.output,
        forward_file = rules.pre_process1.output
    output:
        output_log = config["inputdir"] + "{fq_file}/output/processed/trimmomatic.log",
        output_manifest = config["inputdir"] + "{fq_file}/input/pe-33v1-manifest.tsv"
    conda:
        "environments/multiqc.yaml"
    params:
        forward_output_P = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_1_filtered_P.fastq.gz",
        reverse_output_P = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_2_filtered_P.fastq.gz",
        forward_output_U = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_1_filtered_U.fastq.gz",
        reverse_output_U = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_2_filtered_U.fastq.gz",
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


rule post_fastqc_forward:
    input:
        rules.trimmomatic.output.output_manifest
    output:
        forward_fastqc = config["inputdir"] + "output/post_fastqc/{fq_file}_1_filtered_P_fastqc.html",
        forward_fastqc_zip = config["inputdir"] + "output/post_fastqc/{fq_file}_1_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["inputdir"] + "output/post_fastqc/",
        input_file=config["inputdir"] + "{fq_file}/input/{fq_file}_1.fastq.gz",
        input_file2=rules.trimmomatic.params.forward_output_P
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
rule post_fastqc_reverse:
    input:
        input_file2=rules.post_fastqc_forward.output.forward_fastqc_zip
    output:
        reverse_fastqc=config["inputdir"] + "output/post_fastqc/{fq_file}_2_filtered_P_fastqc.html",
        reverse_fastqc_zip=config["inputdir"] + "output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip"
    conda:
        "environments/multiqc.yaml"
    params:
        output_dir=config["inputdir"] + "output/post_fastqc/",
        input_file=config["inputdir"] + "{fq_file}/input/{fq_file}_2.fastq.gz",
        input_file2=rules.trimmomatic.params.reverse_output_P
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

rule post_multiqc:
    input:
        expand(config["inputdir"] + "output/post_fastqc/{fq_file}_2_filtered_P_fastqc.zip",fq_file=fq_files)
    output:
        config["inputdir"] + "output/post_multiqc/multiqc_report.html"
    conda:
        "environments/multiqc.yaml"
    params:
        outputdir=config["inputdir"] + "output/post_multiqc/",
        inputtdir=config["inputdir"] + "output/post_fastqc/"
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


rule make_cache:
    input:
        rules.post_multiqc.output
    output:
        cache_dir = directory(config["inputdir"] + "temp/{fq_file}"),
        done = config["inputdir"] + "{fq_file}/output/steps/{fq_file}_cache_made.txt"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        """
        mkdir -p {output.cache_dir}
        sbatch bash_scripts/make_cache.sh {output.cache_dir}/cache {output.done}
        python3 {config[tooldir]}wait_file.py {output.done} --seconds=40
        """

rule make_artefact:
    input:
        cache = rules.make_cache.output,
        manifest = rules.trimmomatic.output.output_manifest
    output:
        config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_reads.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    message:
        "@#"
        "Qiime import:   "
        "qiime tools import   "
        "   --type 'SampleData[PairedEndSequencesWithQuality]'   "
        "   --input-path {input.manifest}   "
        "   --input-format PairedEndFastqManifestPhred33   "
        "   --output-path {output}"
    shell:
        'sbatch bash_scripts/create_artifact.sh {input.manifest} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule make_merged:
    input:
        input1 = rules.trimmomatic.output
    output:
        touch(config["inputdir"] + "output/steps/{fq_file}_merged_fw_done.txt")
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        forward_output_P = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_1_filtered_P.fastq.gz",
        reverse_output_P = config["inputdir"] + "{fq_file}/output/processed/{fq_file}_2_filtered_P.fastq.gz",
        fw_output_merged = config["inputdir"] + "output/processed/merged_1.fastq.gz",
        rv_output_merged = config["inputdir"] + "output/processed/merged_2.fastq.gz",
        output_dir = config["inputdir"] + "output/processed"
    message:
        "@#"
        "Merging:   "
        "cat {params.forward_output_P} {params.fw_output_merged}    "
        "cat {params.reverse_output_P} {params.rv_output_merged}    "
        "@#"
    shell:
        'mkdir -p {params.output_dir};'
        'sbatch bash_scripts/merging.sh {params.forward_output_P} {params.fw_output_merged} {params.reverse_output_P} {params.rv_output_merged};'
        'python3 {config[tooldir]}wait_file.py {params.fw_output_merged} {params.rv_output_merged} --seconds=40;'

# i am testing if i can add the ids of the samples to the merged manifest and see if the index will be able to accept it. (if it doesnt change this part back.
# rule make_artefact_merged:
#     input:
#         expand(config["inputdir"] + "output/steps/{fq_file}_merged_fw_done.txt",fq_file=fq_files)
#     output:
#         config["inputdir"] + "output/processed/" + config["naming_convention"] + "merged_reads.qza"
#     conda:
#         "environments/qiime2_metagenome-2024.10.yaml"
#     params:
#         manifest = config["inputdir"] + "output/processed/manifest.txt"
#     shell:
#         'printf "sample-id,absolute-filepath,direction\nCombined_samples,{rules.make_merged.params.fw_output_merged},forward\nCombined_samples,{rules.make_merged.params.rv_output_merged},reverse\n" > {params.manifest};'
#         'sbatch bash_scripts/create_artifact.sh {params.manifest} {output};'
#         'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule make_artefact_merged:
    input:
        expand(config["inputdir"] + "output/steps/{fq_file}_merged_fw_done.txt",fq_file=fq_files)
    output:
        config["inputdir"] + "output/processed/" + config["naming_convention"] + "merged_reads.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        manifest = config["inputdir"] + "output/processed/manifest.txt"
    message:
        "@#"
        "Create manifest and artifact:   "
        "qiime tools import "
        "   --type 'SampleData[PairedEndSequencesWithQuality]' "
        "   --input-path {params.manifest} "
        "   --input-format PairedEndFastqManifestPhred33    "
        "   --output-path {output}    "
        "@#"
    shell:
        'bash bash_scripts/get_manifest_1_0.sh {config[inputdir]} {rules.make_merged.params.fw_output_merged} {rules.make_merged.params.rv_output_merged} {params.manifest};'
        'sbatch bash_scripts/create_artifact.sh {params.manifest} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule megahit:
    input:
        rules.make_artefact_merged.output
    output:
        config["inputdir"] + "output/processed/" + config["naming_convention"] + "contigs.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        #max 10
        thread = 8
    message:
        """@#"
        "Megahit assemble contigs:   "
        "qiime assembly assemble-megahit    "
        "   --i-seqs {input}    "
        "   --p-presets "meta-sensitive"    "
        "   --p-num-cpu-threads {params.thread} "
        "   --o-contigs {output}    "
        "   --verbose"
        "@#"""
    shell:
        'sbatch bash_scripts/megahit.sh {input} {output} {params.thread};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule index_contigs:
    input:
        rules.megahit.output
    output:
        config["inputdir"] + "output/processed/" + config["naming_convention"] + "contigs_indexed.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        #max 8
        thread = 8
    message:
        "@#"
        "Index contigs:   "
        "qiime assembly index-contigs   "
        "   --i-contigs {input} "
        "   --o-index {output} "
        "   --p-threads {params.thread} "
        "-  -verbose"
        "@#"
    shell:
        'sbatch bash_scripts/index_megahit.sh {input} {output} {params.thread};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'


#change
rule mapping_contigs:
    input:
        input_index = rules.index_contigs.output,
        input_sample = rules.make_artefact.output
    output:
         config["inputdir"] + "{fq_file}/output/processed/{fq_file}_" + config["naming_convention"] + "read_aln_contigs.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        #max 8
        thread = 8
    message:
        """@#"
        "Map reads to contigs:   "
        "qiime assembly _map-reads-to-contigs "
        "   --i-index {input.input_index} "
        "   --i-reads  {input.input_sample} "
        "   --o-alignment-map {output} "
        "   --p-threads {params.thread} "
        "   --p-sensitivity "sensitive" 
            --p-seed 711 
            --verbose"
        "@#"""
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {input.input_sample} {output} {params.thread};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule export_data:
    input:
        rules.megahit.output
    output:
        config["inputdir"] + "output/processed/emapper_input/{fq_file}_contigs.fa"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        output_dir = config["inputdir"] + "output/processed/emapper_input"
    message:
        "@#"
        "Exporting megahit contigs: "
        "qiime tools export "
        "   --input-path {input} "
        "   --output-path {params.output_dir} "
        "@#"
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/export.sh {input} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule emapper:
    input:
        rules.export_data.output,
    output:
        fasta = config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna.emapper.genepred.fasta",
        ggf = config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna.emapper.genepred.gff",
        hits =  config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna.emapper.hits"
    conda:
        "environments/eggnog.yaml"
    params:
        #max 8
        thread = 8,
        output_dir = config["inputdir"] + "{fq_file}/output/processed/emapper_output",
        temp_dir = "/export/jippe/temp/{fq_file}/",
        orthologs = config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna.emapper.seed_orthologs",
        moved_dir = config["inputdir"] + "{fq_file}/output/processed/emapper_output/orthologs"
    message:
        """@#
        "eggnog mapping:   "
        "python3 eggnog-mapper/emapper.py "
        "   -i  {input} "
        "   -o merged_dna "
        "   -m diamond "
        "   --dmnd_db "/export/databases/cache/data/979a93e7-e318-437a-9e35-37c981f7d787/data/ref_db.dmnd" "
        "   --itype metagenome "
        "   --output_dir {params.output_dir} "
        "   --temp_dir {params.temp_dir} "
        "   --cpu {params.thread}"
        "@#"""
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'mkdir -p {params.temp_dir};'
        'sbatch bash_scripts/emapper_manual.sh {input} merged_dna {params.output_dir} {params.temp_dir} {params.thread};'
        'python3 {config[tooldir]}wait_file.py {output.fasta} {output.ggf} {output.hits} --seconds=40;'
        'mkdir -p {params.moved_dir};'
        'mv {params.orthologs} {params.moved_dir};'

rule import_hits:
    input:
        # rules.emapper.output.hits
        rules.emapper.output.hits
    output:
        config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna_emapper_hits.qza"
    conda:
        "environments/qiime2_2024.5_moshpit.yaml"
    message:
        "@#"
        "Importing egnog mapper hits:   "
        "qiime tools import "
        "   --input-path {rules.emapper.params.moved_dir} "
        "   --output-path {output} "
        "   --type '"'SampleData[BLAST6]'"'"
        "@#"
    shell:
        'sbatch bash_scripts/import_hits.sh {rules.emapper.params.moved_dir} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule moshpit_eggnog_annotate:
    input:
        rules.import_hits.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna_annotations_contig.qza"
    conda:
        "environments/qiime2_2024.5_moshpit.yaml"
    params:
        eggnog_annot="/export/databases/cache:eggnog_annotation",
        #max = 8
        threads = 8
    message:
        "@#"
        "Annoting the hits:   "
        "qiime moshpit eggnog-annotate "
        "   --i-eggnog-hits {input} "
        "   --i-eggnog-db {params.eggnog_annot} "
        "   --p-num-cpus {params.threads} "
        "   --o-ortholog-annotations {output}"
        "   --verbose"
        "@#"
    shell:
        'sbatch bash_scripts/emapper_annotate.sh {input} {params.eggnog_annot} {params.threads} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule export_data_annotated:
    input:
        rules.moshpit_eggnog_annotate.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/emapper_output/merged_dna.emapper.annotations"
    conda:
        "environments/qiime2_2024.5_moshpit.yaml"
    params:
        output_dir = config["inputdir"] + "{fq_file}/output/processed/emapper_output"
    message:
        "@#"
        "Export mosphit annotated:   "
        "qiime tools export "
        "   --input-path {input} "
        "   --output-path {params.output_dir}"
        "@#"
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/export.sh {input} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule export_to_bam:
    input:
        # signal = rules.export_data_annotated.output,
        inputfile = rules.mapping_contigs.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/alignment/{fq_file}_alignment.bam"
    conda:
        "environments/qiime2_2024.5_moshpit.yaml"
    params:
        output_dir = config["inputdir"] + "{fq_file}/output/processed/alignment/"
    message:
        "@#"
        "Export to bam:   "
        "qiime tools export "
        "   --input-path {input.inputfile} "
        "   --output-path {params.output_dir}"
        "@#"
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/export.sh {input.inputfile} {params.output_dir};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'


rule sort_bam:
    input:
        rules.export_to_bam.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/alignment/{fq_file}_sorted_alignment.bam"
    conda:
        "environments/qiime2_2024.5_moshpit.yaml"
    message:
        "@#"
        "sort to bam: "
        "samtools sort "
        "   -n "
        "   -o {output} "
        "   {input}"
        "@#"
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/sort_bam.sh {input} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

# rule sort_bam_index:
#     input:
#         rules.sort_bam.output
#     output:
#         config["inputdir"] + "{fq_file}/output/processed/alignment/{fq_file}_sorted_alignment.bam.bai"
#     conda:
#         "environments/qiime2_2024.5_moshpit.yaml"
#     params:
#         cpu = 12
#     message:
#         "@#"
#         "Index bam:   "
#         "samtools index "
#         "   -@ {params.cpu} "
#         "   -o {output} "
#         "   {input}"
#         "@#"
#     shell:
#         # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
#         'sbatch bash_scripts/bam_index.sh {params.cpu} {output} {input};'
#         'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule htseq_count:
    input:
        aln_bam = rules.sort_bam.output,
        gff_file = rules.emapper.output.ggf
    output:
        config["inputdir"] + "{fq_file}/output/processed/{fq_file}_merged_HTSeq_gene_counts.txt"
    conda:
        "environments/multiqc.yaml"
    message:
        "@#"
        "   "
        "htseq-count "
        "   -f bam "
        "   -r name "
        "   -s no "
        "   -t CDS "
        "   -i ID {input.aln_bam} {input.gff_file} > {output} "
        "@#"
    shell:
        # 'sbatch bash_scripts/map_reads_to_contigs.sh {input.input_index} {rules.trimmomatic.params.forward_output_P} {rules.trimmomatic.params.reverse_output_P} {output} {params.thread};'
        'sbatch bash_scripts/htseq_count.sh  {input.aln_bam} {input.gff_file} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'
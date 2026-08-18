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
    print("\n✔ Detected sample/input/*.fastq.gz → skipping renaming step")
else:
    print("\n🔄 Organizing FASTQ files into sample/input/")
    change_names(folder)

# 📂 List samples
fq_files = list_sample_folders(folder)
print("\nDetected sample folders:")
for s in fq_files:
    print(f" - {s}")



conf_thresh = str(config["CONF_THRESH"] )

rule all:
    input:
        expand(config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_table.qzv", conf_thresh=conf_thresh),
        expand(config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_barplot.qzv", conf_thresh=conf_thresh),
        expand(config["inputdir"] + "output/artefacts/report_diamond_pathogens_plant/" + config["naming_convention"] + "{fq_file}_diaomond_pathogens_plant.f6", fq_file=fq_files),
        expand(config["inputdir"] + "output/artefacts/report_diamond_pathogens_human/" + config["naming_convention"] + "{fq_file}_diaomond_pathogens_human.f6", fq_file=fq_files)
# # Old Steps
# rule Sorting:
#     output:
#         forward_file = config["inputdir"] + "raw_data/{fq_file}_1.fq.gz",
#         reverse_file = config["inputdir"] + "raw_data/{fq_file}_2.fq.gz"
#     conda:
#         "environments/multiqc.yaml"
#     params:
#         forward_file =  config["inputdir"] + "{fq_file}_1.fq.gz",
#         reverse_file = config["inputdir"] + "{fq_file}_2.fq.gz"
#     message:
#         "@#"
#         "move to raw_data folder:   "
#         "mv "
#         "   *_1.fq.gz raw_data/*_1/2.fq.gz"
#         "@#"
#     shell:
#        'mv {params.forward_file} {output.forward_file}; mv {params.reverse_file} {output.reverse_file}'


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

rule combine_fw_rv:
    input:
        rules.trimmomatic.output
    output:
        config["inputdir"] + "{fq_file}/output/processed/{fq_file}_merged.fastq.gz"
    conda:
        "environments/multiqc.yaml"
    params:
        threads=1,
        fw = rules.trimmomatic.params.forward_output_P,
        rv = rules.trimmomatic.params.reverse_output_P
    message:
        "@#"
        "Quality control forward:   "
        "fastqc "
        "@#"
    shell:
        "cat {params.fw} {params.rv} > {output}"

rule pathogen_human:
    input:
        rules.combine_fw_rv.output
    output:
        config["inputdir"] + "output/artefacts/report_diamond_pathogens_human/" + config["naming_convention"] + "{fq_file}_diaomond_pathogens_human.f6"
    conda:
        "environments/multiqc.yaml"
    params:
        threads=16,
        db="/export/databases/Wetsus/MetaPatho/human/MetaPatho-HBPDB.v1.0"
    message:
        "@#"
        "diamond blastx --db $1 --query $2 --outfmt 6 --threads $3 --max-target-seqs 1 --quiet -e 1e-10 --more-sensitive --block-size 10 --query-cover 70 --out $4"
        "@#"
    shell:
        "sbatch bash_scripts/pathogen.sh {params.db} {input} {params.threads} {output};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=40"

rule pathogen_plant:
    input:
        rules.combine_fw_rv.output
    output:
        config["inputdir"] + "output/artefacts/report_diamond_pathogens_plant/" + config["naming_convention"] + "{fq_file}_diaomond_pathogens_plant.f6"
    conda:
        "environments/multiqc.yaml"
    params:
        threads=16,
        db="/export/databases/Wetsus/MetaPatho/plant/MetaPatho-PBPDB.v1.0"
    message:
        "@#"
        "diamond blastx --db $1 --query $2 --outfmt 6 --threads $3 --max-target-seqs 1 --quiet -e 1e-10 --more-sensitive --block-size 10 --query-cover 70 --out $4"
        "@#"
    shell:
        "sbatch bash_scripts/pathogen.sh {params.db} {input} {params.threads} {output};"
        "python3 {config[tooldir]}wait_file.py {output} --seconds=40"


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
    shell:
        'sbatch bash_scripts/create_artifact.sh {input.manifest} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule classify_kraken:
    input:
        artefact = rules.make_artefact.output,
        cache_dir = rules.make_cache.output.cache_dir
    output:
        report = config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_report.qza",
        hits = config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_hit.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    params:
        input_file = lambda wildcards, input: f"{input.cache_dir}/cache/data/{config['naming_convention']}{wildcards.fq_file}_reads.qza"
    shell:
        'cp -au {input.artefact} {input.cache_dir}/cache/data/;'
        'sbatch bash_scripts/classify_kraken.sh {params.input_file} {config[kraken2_db]} {config[CONF_THRESH]} {output.report} {output.hits};'
        'python3 {config[tooldir]}wait_file.py {output.report} {output.hits} --seconds=40;'


rule kraken_features:
    input:
        report = rules.classify_kraken.output.report,
        hits = rules.classify_kraken.output.hits
    output:
        table = config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_table.qza",
        taxonomy = config["inputdir"] + "{fq_file}/output/artefacts/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_taxonomy.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/kraken_features.sh {input.report} {output.table} {output.taxonomy};'
        'python3 {config[tooldir]}wait_file.py {output.table} {output.taxonomy} --seconds=40;'


rule pre_collate_features:
    input:
        order = rules.kraken_features.output,
        input_file = rules.classify_kraken.output.report
    output:
        config["inputdir"] + "output/artefacts/reports/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_report.qza"
    shell:
        "cp {input.input_file} {output}"


rule collate_features:
    input:
        expand(config["inputdir"] + "output/artefacts/reports/" + config["naming_convention"] + "{fq_file}_{conf_thresh}_kraken_report.qza", fq_file=fq_files, conf_thresh=conf_thresh),
    output:
        config["inputdir"] + "output/artefacts/" + config["naming_convention"] + "{conf_thresh}_collated_kraken_reports.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/collate_features.sh {input} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule estimate_bracken:
    input:
        rules.collate_features.output
    output:
        report = config["inputdir"] + "output/artefacts/bracken_reports/qza/" + config["naming_convention"] + "{conf_thresh}_bracken_species_report.qza",
        taxonomy = config["inputdir"] + "output/artefacts/bracken_reports/qza/" + config["naming_convention"] + "{conf_thresh}_bracken_species_taxonomy.qza",
        table = config["inputdir"] + "output/artefacts/bracken_reports/qza/" + config["naming_convention"] + "{conf_thresh}_bracken_species_table.qza"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/estimate_bracken.sh {input} {config[bracken_db]} {config[read_length]} {output.report} {output.taxonomy} {output.table};'
        'python3 {config[tooldir]}wait_file.py {output.report} {output.taxonomy} {output.table} --seconds=40;'

rule summarize_table:
    input:
        rules.estimate_bracken.output.table
    output:
        config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_table.qzv"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/summarize_table.sh {input} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'

rule export_bracken_taxonomy:
    input:
        report = rules.estimate_bracken.output.report,
        taxonomy = rules.estimate_bracken.output.taxonomy
    output:
        report = directory(config["inputdir"] + "output/artefacts/bracken_reports/" + config["naming_convention"] + "{conf_thresh}_bracken_species_report_output"),
        taxonomy = directory(config["inputdir"] + "output/artefacts/bracken_reports/" + config["naming_convention"] + "{conf_thresh}_bracken_species_table_output")
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/export_tool.sh {input.report} {output.report};'
        'sbatch bash_scripts/export_tool.sh {input.taxonomy} {output.taxonomy};'
        'python3 {config[tooldir]}wait_file.py {output.report} {output.taxonomy} --seconds=40;'

rule taxa_barplot:
    input:
        order1 = rules.export_bracken_taxonomy.output.report,
        order2 = rules.export_bracken_taxonomy.output.taxonomy,
        table = rules.estimate_bracken.output.table,
        taxonomy = rules.estimate_bracken.output.taxonomy
    output:
        config["inputdir"] + "output/artefacts/bracken_reports/qzv/" + config["naming_convention"] + "{conf_thresh}_bracken_species_barplot.qzv"
    conda:
        "environments/qiime2_metagenome-2024.10.yaml"
    shell:
        'sbatch bash_scripts/taxa_barplot.sh {input.table} {input.taxonomy} {output};'
        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'


# rule bbmap_interleave:
#     input:
#         input1 = rules.trimmomatic.output.output_log,
#         input2 = rules.post_multiqc.output
#     output:
#         startdir + "{dataset}/output/processed/interleaved.fastq.gz"
#     conda:
#         "multiqc.yaml"
#     params:
#         forward_output = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz",
#         reverse_output = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz",
#     message:
#         "@#"
#         "Interleaving:   "
#         "reformat.sh"
#         "   in1=*/output/processed/forward_filtered_P.fastq.gz"
#         "   in2=*/output/processed/reverse_filtered_P.fastq.gz"
#         "   out=*/output/processed/interleaved.fastq.gz"
#         "@#"
#     shell:
#        'sbatch bash_scripts/bbmap_interleave.sh {params.forward_output} {params.reverse_output} {output};'
#        'python3 {config[tooldir]}wait_file.py {output} --seconds=40;'
#        'rm {params.forward_output}; rm {params.reverse_output}'
#
# rule classify_kma:
#     input:
#         rules.bbmap_interleave.output
#     output:
#         output1 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.fsa",
#         output2 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.frag.gz",
#         output3 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.aln",
#         output4 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.res",
#         output5 = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output.tsv"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_name = startdir + "{dataset}/output/classified/card/kma_card_3.2.7_output"
#     message:
#         "@#"
#         "Classifying card:   "
#         "kma "
#         "   -int $1"
#         "   -o $2"
#         "   -t_db /export/databases/Wetsus/Card/3.2.7/kma_index/nucleotide_fasta_protein_homolog_model_variants"
#         "   -tsv"
#         "   -1t1"
#         "   -cge"
#         "   -hmm"
#         "   -ex_mode"
#         "   -mem_mode"
#         "@#"
#     shell:
#         'sbatch bash_scripts/kma_classify.sh {input} {params.output_name};'
#         'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
#
# rule classify_kma_resfinder:
#     input:
#         rules.bbmap_interleave.output
#     output:
#         output1 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.fsa",
#         output2 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.frag.gz",
#         output3 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.aln",
#         output4 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.res",
#         output5 = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output.tsv"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_name = startdir + "{dataset}/output/classified/resfinder/kma_resfinder_3.2.7_output"
#     message:
#         "@#"
#         "Classifying Resfinder:   "
#         "kma "
#         "   -int $1"
#         "   -o $2"
#         "   -t_db /export/databases/Wetsus/Resfinder2/ResFinder2"
#         "   -tsv"
#         "   -1t1"
#         "   -cge"
#         "   -hmm"
#         "   -ex_mode"
#         "   -mem_mode"
#         "@#"
#     shell:
#         'sbatch bash_scripts/kma_classify_resfinder.sh {input} {params.output_name};'
#         'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
#
# rule classify_kma_silva:
#     input:
#         rules.bbmap_interleave.output
#     output:
#         output1 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.fsa",
#         output2 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.frag.gz",
#         output3 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.aln",
#         output4 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.res",
#         output5 = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output.tsv"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_name = startdir + "{dataset}/output/classified/silva/kma_silva_138.2_output"
#     message:
#         "@#"
#         "Classifying Resfinder:   "
#         "kma "
#         "   -int $1"
#         "   -o $2"
#         "   -t_db /export/databases/Silva/138.2/kma/silva_138.2_nr99_kma"
#         "   -tsv"
#         "   -1t1"
#         "   -cge"
#         "   -hmm"
#         "   -ex_mode"
#         "   -mem_mode"
#         "@#"
#     shell:
#         'sbatch bash_scripts/kma_classify_silva138.2.sh {input} {params.output_name};'
#         'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
#
# rule classify_kma_panres:
#     input:
#         rules.bbmap_interleave.output
#     output:
#         output1=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.fsa",
#         output2=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.frag.gz",
#         output3=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.aln",
#         output4=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.res",
#         output5=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output.tsv"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_name=startdir + "{dataset}/output/classified/panres/panres_31_07_2025_output"
#     message:
#         "@#"
#     "Classifying Resfinder:   "
#     "kma "
#     "   -int $1"
#     "   -o $2"
#     "   -t_db /export/databases/Wetsus/Panres/kma/panres_31_07_2025"
#     "   -tsv"
#     "   -1t1"
#     "   -cge"
#     "   -hmm"
#     "   -ex_mode"
#     "   -mem_mode"
#     "@#"
#     shell:
#         'sbatch bash_scripts/kma_panres.sh {input} {params.output_name};'
#         'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
#
#
# rule classify_kma_pr2:
#     input:
#         rules.bbmap_interleave.output
#     output:
#         output1=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.fsa",
#         output2=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.frag.gz",
#         output3=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.aln",
#         output4=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.res",
#         output5=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output.tsv"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_name=startdir + "{dataset}/output/classified/pr2/pr2_5.1.0_output"
#     message:
#         "@#"
#     "Classifying Resfinder:   "
#     "kma "
#     "   -int $1"
#     "   -o $2"
#     "   -t_db /export/databases/Wetsus/PR2/5.1.0/kma_16/pr2_5.1.0_SSU"
#     "   -tsv"
#     "   -1t1"
#     "   -cge"
#     "   -hmm"
#     "   -ex_mode"
#     "   -mem_mode"
#     "@#"
#     shell:
#         'sbatch bash_scripts/kma_pr2.sh {input} {params.output_name};'
#         'python3 {config[tooldir]}wait_file.py {output.output1} {output.output2} {output.output3} {output.output4} --seconds=1200;'
#
#



























##############################
#
# rule trimmomatic:
#     input:
#         rules.pre_process2.output
#     output:
#         trimmomatic_output
#         output_log = startdir + "{dataset}/output/processed/trimmomatic.log",
#         output_2 = startdir + "{dataset}/output/done/3_trimmomatic_done.txt"
#     conda:
#         "multiqc.yaml"
#     params:
#         forward_file=startdir + "{dataset}/output/processed/forward_filtered.fastq.gz",
#         reverse_file=startdir + "{dataset}/output/processed/reverse_filtered.fastq.gz",
#         forward_output = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz",
#         reverse_output = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz",
#         forward_output_U = startdir + "{dataset}/output/processed/forward_filtered_U.fastq.gz",
#         reverse_output_U = startdir + "{dataset}/output/processed/reverse_filtered_U.fastq.gz"
#     message:
#         "@#"
#         "Trimming:   "
#         "java -jar /export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar   "
#         "   PE"
#         "   -threads 1"
#         "   */output/processed/forward_filtered.fastq.gz */output/processed/reverse_filtered.fastq.gz */output/processed/forward_filtered_P.fastq.gz */output/processed/forward_filtered_U.fastq.gz */output/processed/reverse_filtered_P.fastq.gz */output/processed/reverse_filtered_U.fastq.gz"
#         "   ILLUMINACLIP:/export/jippe/jsil/programs/Transcriptoom/versions/1.0/tools/trimmomatic/Trimmomatic-0.39/adapters/Adapters_Macrogen.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:30 HEADCROP:1 2> */output/processed/trimmomatic.log"
#         "@#"
#     shell:
#        'bash bash_scripts/continuing.sh;'
#        'sbatch bash_scripts/trimmomatic.sh {params.forward_file} {params.reverse_file} {startdir}{wildcards.dataset}/output/processed {output.output_log};'
#        'python3 {config[tooldir]}wait_file.py {params.forward_output} {params.reverse_output} --seconds=200;'
#        'rm {params.forward_output_U}; rm {params.reverse_output_U};'
#        'rm {params.forward_file}; rm {params.reverse_file};'
#        'touch {output.output_2}'

# rule post_fastqc_forward:
#     input:
#         trimmomatic = rules.trimmomatic.output.output_log
#     output:
#         forward_fastqc = startdir + "{dataset}/output/post_fastqc/forward_filtered_P_fastqc.html",
#         forward_fastqc_zip = startdir + "{dataset}/output/post_fastqc/forward_filtered_P_fastqc.zip"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_dir = startdir + "{dataset}/output/post_fastqc/",
#         forward_file = startdir + "{dataset}/output/processed/forward_filtered_P.fastq.gz"
#     message:
#         "@#"
#         "Post trimmomatic Quality control forward:   "
#         "fastqc "
#         "   */output/processed/forward_filtered_P.fastq.gz "
#         "   -o */output/post_fastqc/ "
#         "@#"
#     shell:
#        'sbatch bash_scripts/fastq.sh {params.forward_file} {params.output_dir};'
#        'python3 {config[tooldir]}wait_file.py {output.forward_fastqc} {output.forward_fastqc_zip} --seconds=40;'
#
# rule post_fastqc_reverse:
#     input:
#         fw_fastqc = rules.post_fastqc_forward.output.forward_fastqc,
#         fw_fastqc_zip = rules.post_fastqc_forward.output.forward_fastqc_zip
#     output:
#         reverse_fastqc = startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.html",
#         reverse_fastqc_zip = startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.zip"
#     conda:
#         "multiqc.yaml"
#     params:
#         output_dir = startdir + "{dataset}/output/post_fastqc/",
#         reverse_file = startdir + "{dataset}/output/processed/reverse_filtered_P.fastq.gz"
#     message:
#         "@#"
#         "Post trimmomatic Quality control reverse:   "
#         "fastqc "
#         "   */output/processed/reverse_filtered_P.fastq.gz "
#         "   -o */output/post_fastqc/ "
#         "@#"
#     shell:
#        'sbatch bash_scripts/fastq.sh {params.reverse_file} {params.output_dir};'
#        'python3 {config[tooldir]}wait_file.py {output.reverse_fastqc} {output.reverse_fastqc_zip}  --seconds=40;'
#
# rule post_multiqc:
#     input:
#         elly = expand(startdir + "{dataset}/output/post_fastqc/reverse_filtered_P_fastqc.zip", dataset=input_wildcard)
#     output:
#         startdir + "post_multiqc/multiqc_report.html"
#     conda:
#         "multiqc.yaml"
#     params:
#         outputdir = startdir + "post_multiqc/"
#     message:
#         "@#"
#         "Post processing Multiqc:   "
#         "multiqc "
#         "   */output/fastqc/* "
#         "   --outdir post_multiqc/ "
#         "   -f "
#         "   -d "
#         "   -s"
#         "@#"
#     shell:
#          "multiqc {startdir}*/output/post_fastqc/* --outdir {params.outputdir} -f -d -s"
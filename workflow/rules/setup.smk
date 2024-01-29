from os import listdir, path
from snakemake.utils import min_version


# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

V2_MURATES_URL = config["urls"]["v2_murates"]
V3_CONSTRAINT_URL = config["urls"]["v3_constraint"]
V2_EXOMES = config["urls"]["v2_exomes"]
V3_GENOMES = config["urls"]["v3_genomes"]
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]

# ------------- #
# I/O           #
# ------------- #

# gnomad v2 supplement data
GNOMAD_V2_SUPPLEMENT = path.join(
    INSTALL_DIR, "supplement", "41586_2020_2308_MOESM4_ESM.zip"
)
GNOMAD_V2_MURATES = path.join(
    INSTALL_DIR, "supplement", "gnomad_v2.supplement-f10.murates.tsv"
)

# gnomadV3 nocoding constraint data
GNOMAD_V3_SUPPLEMENT = path.join(INSTALL_DIR, "supplement", "media-2.zip?download=true")
GNOMAD_V3_CONSTRAINT = path.join(
    INSTALL_DIR, "supplement", "gnomad_v3.supplement-f03.constraintZ.100bp.hg38.tsv"
)


# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        GNOMAD_V2_MURATES,
        GNOMAD_V3_CONSTRAINT,


rule download_gnomad_v2_supplement:
    message:
        """
        Downloads gnomAD v2 supplement.
        """
    output:
        temp(GNOMAD_V2_SUPPLEMENT),
    params:
        url=V2_MURATES_URL,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/download_gnomad_supplement.stdout",
        stderr="workflow/logs/download_gnomad_supplement.stderr",
    shell:
        "wget {params.url} -O {output}"


rule extract_gnomad_v2_murates:
    message:
        """
        Unpacks supllement and extracts mutation rates table.
        """
    input:
        rules.download_gnomad_v2_supplement.output,
    output:
        GNOMAD_V2_MURATES,
    params:
        target="supplement/supplementary_dataset_10_mutation_rates.tsv.gz",
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/extract_v2_murates.stdout",
        stderr="workflow/logs/extract_v2_murates.stderr",
    shell:
        "unzip -p {input} {params.target} | zcat > {output}"


rule download_gnomad_v3_noncoding_constraint_supplement:
    message:
        """
        Downloads gnomAD v3 noncoding constraint supplement.
        """
    output:
        temp(GNOMAD_V3_SUPPLEMENT),
    params:
        url=V3_CONSTRAINT_URL,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/download_gnomad_v3_noncoding_constraint_supplement.stdout",
        stderr="workflow/logs/download_gnomad_v3_noncoding_constraint_supplement.stderr",
    shell:
        "wget {params.url} -O {output}"


rule extract_gnomad_v3_noncoding_constraint:
    message:
        """
        Extracts noncoding constraint data from supplement.
        """
    input:
        rules.download_gnomad_v3_noncoding_constraint_supplement.output,
    output:
        GNOMAD_V3_CONSTRAINT,
    conda:
        "../envs/gnomad-maps.yaml"
    threads: 1
    log:
        stdout="workflow/logs/extract_v2_murates.stdout",
        stderr="workflow/logs/extract_v2_murates.stderr",
    shell:
        "unzip -p {input} Supplementary_Datasets/Supplementary_Data_3.bed.gz | zcat > {output}"


# ---------------------------------------------------- #
# Downloading genetic data - leave for example         #
# ---------------------------------------------------- #
# rule download_gnomad_v2_exome:
#     message:
#         """
#         VCF sites download
#         """
#     output:
#         dir("resources/data/gnomad/variants/v2/exome"),
#     params:
#         url=V2_EXOMES,
#     conda:
#         "gnomad"
#     threads: 1
#     log:
#         stdout="workflow/logs/download_gnomad_v2_exome.stdout",
#         stderr="workflow/logs/download_gnomad_v2_exome.stderr",
#     shell:
#         """
#         gsutil cp {params.url}.tbi {output} &&
#         gsutil cp {params.url}.bgz {output}
#         """
# rule download_gnomad_v3_genome:
#     message:
#         """
#         Hail table download
#         """
#     output:
#         dir("resources/data/gnomad/variants/v3/genome"),
#     params:
#         url=V3_GENOMES,
#     conda:
#         "gnomad"
#     threads: 1
#     log:
#         stdout="workflow/logs/download_gnomad_v3_genome.stdout",
#         stderr="workflow/logs/download_gnomad_v3_genome.stderr",
#     shell:
#         """
#         gsutil cp -r {params.url}i {output}
#         """

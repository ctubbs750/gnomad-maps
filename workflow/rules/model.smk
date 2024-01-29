from os import listdir, path
from snakemake.utils import min_version


# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

VARIANTS = config["variants"]
CLASSES = config["classes"]
GENOME = config["genome"]
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]

# ------------- #
# I/O           #
# ------------- #

# gnomad v2 supplement data
CLASS_SNVS = path.join(PROCESS_DIR, "snvs", "gnomad_v3.snvs.{variant_class}.bed")
TRICONTEXT = path.join(
    PROCESS_DIR, "snvs", "gnomad_v3.snvs.{variant_class}.tricontext.bed"
)
MURATES = path.join(
    PROCESS_DIR, "snvs", "gnomad_v3.snvs.{variant_class}.tricontext.murates.bed"
)

# Calibrated maps
CALIBRATED_MAPS = path.join(PROCESS_DIR, "maps", "maps-calibrated.pickle")

# ------------- #
# Rules         #
# ------------- #


wildcard_constraints:
    variant_class="\w+",


rule all:
    input:
        expand(MURATES, variant_class=CLASSES),
        CALIBRATED_MAPS,


rule class_snvs:
    message:
        """
        Combs through all gnomAD variants and pulls out variants with VEP of a given class.
        """
    input:
        VARIANTS,
    output:
        temp(CLASS_SNVS),
    params:
        vc=lambda wc: wc.variant_class,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/class_snvs-{variant_class}.stdout",
        stderr="workflow/logs/class_snvs-{variant_class}.stderr",
    shell:
        """
        set +o pipefail
        zcat {input} | vawk 'index($9, "{params.vc}")' > {output}
        """


rule annotate_tricontext_class_snvs:
    message:
        """
        Annotates class snvs with tricontext 
        """
    input:
        rules.class_snvs.output,
    output:
        temp(TRICONTEXT),
    params:
        vc=lambda wc: wc.variant_class,
        genome=GENOME,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/annotate_tricontext_class_snvs-{variant_class}.stdout",
        stderr="workflow/logs/annotate_tricontext_class_snvs-{variant_class}.stderr",
    script:
        "../scripts/tricontext.py"


rule annotate_murates_class_snvs:
    message:
        """
        Annotates class snvs with mutation rates. 
        """
    input:
        rules.annotate_tricontext_class_snvs.output,
        murates="resources/data/gnomad/supplement/gnomad_v2.supplement-f10.murates.tsv",
    output:
        MURATES,
    params:
        vc=lambda wc: wc.variant_class,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/annotate_murates_class_snvs-{variant_class}.stdout",
        stderr="workflow/logs/annotate_murates_class_snvs-{variant_class}.stderr",
    script:
        "../scripts/murates.py"


rule calibrate_maps:
    message:
        """
        Calibrate maps on synvars 
        """
    input:
        MURATES.format(variant_class="synonymous_variant"),
    output:
        CALIBRATED_MAPS,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/calibrate_maps-synonymous_variant.stdout",
        stderr="workflow/logs/calibrate_maps-synonymous_variant.stderr",
    script:
        "../scripts/calibrate.py"

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

# Calculated maps
VARIANT_MAPS = path.join(PROCESS_DIR, "maps", "maps-{variant_class}.tsv")

# Combined maps
COMBINED_MAPS = path.join(PROCESS_DIR, "maps", "maps-vep_classes.tsv")

# ------------- #
# Rules         #
# ------------- #


wildcard_constraints:
    variant_class="\w+",


rule all:
    input:
        COMBINED_MAPS,


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
        temp(MURATES),
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


rule calculate_class_maps:
    message:
        """
        Calibrate maps on all VEP vars 
        """
    input:
        variants=rules.annotate_murates_class_snvs.output,
        model=rules.calibrate_maps.output,
    output:
        temp(VARIANT_MAPS),
    params:
        variant_class=lambda wc: wc.variant_class,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/calculate_class_maps-{variant_class}.stdout",
        stderr="workflow/logs/calculate_class_maps-{variant_class}.stderr",
    script:
        "../scripts/maps.py"


rule combine_class_maps:
    message:
        """
        Combines maps score into single matrix
        """
    input:
        expand(rules.calculate_class_maps.output, variant_class=CLASSES),
    output:
        COMBINED_MAPS,
    conda:
        "../envs/gnomad-maps.yaml"
    log:
        stdout="workflow/logs/combine_class_maps.stdout",
        stderr="workflow/logs/combine_class_maps.stderr",
    shell:
        "head -n 1 {input[0]} > {output}; tail -n +2 -q {input} >> {output}"

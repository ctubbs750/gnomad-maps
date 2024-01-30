"""d"""

from pandas import DataFrame, read_csv
from statsmodels.formula.api import wls

# Snakemake
SYNONYMOUS_VARIANTS = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def neutral_model(matrix: DataFrame):
    """Neutral model, calibrate on synonomous variants"""
    # Return fitted model
    variant_count = list(matrix.constant)
    return wls("sp ~ mu_snp", weights=variant_count, data=matrix).fit()


def main() -> None:
    """Main"""

    # Setup for fast read in
    fields = ["ac", "vep", "context", "mu_snp"]
    dtypes = [int, str, str, float]

    # Read in all variants
    variants = read_csv(
        SYNONYMOUS_VARIANTS,
        sep="\t",
        engine="c",
        usecols=fields,
        dtype=dict(zip(fields, dtypes)),
    )

    # Update with singleton flag
    variants["singleton"] = [1 if ac == 1 else 0 for ac in variants["ac"]]

    # Update with constant
    variants["constant"] = 1

    # Aggregate data to set up for model, count singletons, variants, and mean musnp by factor
    grouped = variants.groupby("vep")
    matrix = grouped.agg({"mu_snp": "mean", "singleton": "sum", "constant": "size"})

    # Combine into model matrix, assign singleton proportion, here "constant" is proxy for variant count
    matrix["sp"] = matrix["singleton"] / matrix["constant"]

    # Calibrate and save model
    neutral_model(matrix=matrix).save(OUTPUT)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()

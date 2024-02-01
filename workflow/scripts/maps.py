"""d"""

from numpy import sqrt, nansum
from pandas import DataFrame, read_csv
from statsmodels.iolib.smpickle import load_pickle

# Snakemake
VARIANTS = snakemake.input[0]  # type: ignore
NEUTRAL_MODEL = snakemake.input[1]  # type: ignore
CLASS = snakemake.params[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def main() -> None:
    """Reads parsed gnomAD variants to DataFrame"""
    # Store data from chunk
    n_singletons = 0
    n_variants = 0
    mu_snp = 0

    # Setup for fast read in
    fields = [6, 12]
    names = ["ac", "mu_snp"]
    dtypes = [int, float]

    # Read chunks of variants
    chunksize = 100000
    for variants in read_csv(
        VARIANTS,
        sep="\t",
        engine="c",
        header=None,
        names=names,
        usecols=fields,
        dtype=dict(zip(fields, dtypes)),
        chunksize=chunksize,
    ):
        # Update with singleton flag
        variants["singleton"] = [1 if ac == 1 else 0 for ac in variants["ac"]]

        # Update with constant
        variants["constant"] = 1

        # Add number of singletons and number of variants (constant) and mu_snp
        n_singletons += sum(variants["singleton"])
        n_variants += len(variants)
        mu_snp += nansum(variants["mu_snp"])

    # Calculate mean mu_snp
    mean_mu_snp = mu_snp / n_variants

    # Create data frame from n_singletons, n_variants, and mu_snp
    matrix = DataFrame(
        {"singleton": [n_singletons], "constant": [n_variants], "mu_snp": [mean_mu_snp]}
    )

    # Update with singleton proportion
    matrix["observed_sp"] = matrix["singleton"] / matrix["constant"]

    # Expected SP
    neutral_model = load_pickle(NEUTRAL_MODEL)

    # Add expected SP
    matrix["expected_sp"] = neutral_model.predict(matrix).item()

    # Adjusted SP
    matrix["maps"] = matrix["observed_sp"] - matrix["expected_sp"]

    # MAPS
    matrix["maps_sem"] = sqrt(
        matrix["observed_sp"] * (1 - matrix["observed_sp"]) / matrix["constant"]
    )

    # Add class
    matrix["class"] = CLASS

    # Save matrix
    matrix.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()

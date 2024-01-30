"""d"""

from pandas import read_csv, merge, DataFrame

# Snakemake
VARIANTS = snakemake.input[0]  # type: ignore
MURATES = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

###
# Functions
###


def read_murates(filepath: str) -> DataFrame:
    """Returns gnomAD murate table as pandas df"""
    fields = ["ref", "alt", "context", "mu_snp"]
    dtypes = [str, str, str, float]
    return read_csv(
        filepath, sep="\t", usecols=fields, engine="c", dtype=dict(zip(fields, dtypes))
    )


def main():
    """Returns variant matrix updates with murates"""
    # Setup murate data
    murate_df = read_murates(MURATES)

    # Params for fast read in
    fields = [
        "chrm",
        "pos0",
        "pos1",
        "ref",
        "alt",
        "af",
        "ac",
        "an",
        "vep",
        "cadd_raw",
        "cadd_phred",
        "context",
    ]
    dtypes = [str, int, int, str, str, float, int, int, str, float, float, str]

    # Read variants
    chunksize = 10000
    for chunk in read_csv(
        VARIANTS,
        sep="\t",
        engine="c",
        header=None,
        names=fields,
        dtype=dict(zip(fields, dtypes)),
        chunksize=chunksize,
    ):
        # Merge datasets
        merge_on = ["ref", "alt", "context"]
        chunk = merge(chunk, murate_df, on=merge_on, how="left", ignore_index=True)

        # Write out
        chunk.to_csv(OUTPUT, sep="\t", index=False, mode="a")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()

"""d"""

from pandas import Series, read_csv
from pysam import FastaFile  # type: ignore

# Snakemake
VARIANTS = snakemake.input[0]  # type: ignore
VEPCLASS = snakemake.params[0]  # type: ignore
GENOME = snakemake.params[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

###
# Functions
###


def genome_fasta(genome_file: str) -> FastaFile:
    """Returns FastaFile class for input genome file"""
    return FastaFile(genome_file)


def genome_sequence(chrm: str, pos0: int, pos1: int, fasta: FastaFile) -> str:
    """Returns string of reference genome sequence for a given interval"""
    # Fetch sequence
    return fasta.fetch(chrm, pos0, pos1).upper()


def fetch_tricontext(variant: Series, genome: FastaFile) -> str:
    """Returns trinucleotide context for input SNV"""
    # Work on copy and set up variant info
    info = variant.copy()
    chrm = info.chrm
    pos0 = info.pos0
    pos1 = info.pos1

    # Set trinucleotide bounds
    l_bound = pos0 - 1
    r_bound = pos1 + 1

    # Query reference sequence for context
    return genome_sequence(chrm, l_bound, r_bound, genome)


def main():
    """Main program"""
    # Setup genome
    genome_fa = genome_fasta(GENOME)

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
    ]
    dtypes = [str, int, int, str, str, float, int, int, str, float, float]

    # Read variants
    chunksize = 10000
    for chunk in read_csv(
        VARIANTS,
        sep="\t",
        header=None,
        engine="c",
        names=fields,
        dtype=dict(zip(fields, dtypes)),
        chunksize=chunksize,
    ):
        # Annotate with tricontext
        chunk["context"] = chunk.apply(
            lambda row: fetch_tricontext(row, genome_fa),
            axis=1,
        )

        # Write out
        chunk.to_csv(OUTPUT, index=False, sep="\t", mode="a", header=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()

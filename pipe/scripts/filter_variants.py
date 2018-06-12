#!/usr/bin/env python3
"""Filter variants based on depth and uniqueness.

To pass this filtering step, variant positions must:
    1. Have min_depth < depth < max_depth across all samples
    2. Show no evidence for being present across multiple samples.
    3. Be a novel heterozygous variant (not loss of heterozygosity)

Variant tables are read into memory directly so only potential variant sites
should be provided to this step.
"""
import sys

import click
from fmbiopy.df import df_subtract
from numpy import float64
from pandas import read_csv


@click.command()
@click.option(
    "--extreme-depth",
    type=click.File("r"),
    help="File containing list of sites with extreme read depth",
)
@click.option(
    "--shared-alleles",
    type=click.File("r"),
    help=(
        "File containing list of potential variants which are found in "
        "multiple samples"
    ),
)
@click.argument("variants", nargs=1, type=click.File("r"))
def filter_variants(extreme_depth, shared_alleles, variants):
    filtered_variants = read_csv(variants, sep="\t")

    for depth_chunk in read_csv(
        extreme_depth,
        delim_whitespace=True,
        chunksize=100000,
        dtype={"CHROM": str, "POS": float64},
    ):
        filtered_variants = df_subtract(filtered_variants, depth_chunk)

    for shared_chunk in read_csv(
        shared_alleles,
        delim_whitespace=True,
        chunksize=100000,
        dtype={"CHROM": str, "POS": float64, "REF": str, "ALT": str},
    ):
        if "snps" in variants:
            snps = shared_chunk[[shared_chunk["ALT"] != "NA"]]
            filtered_variants = df_subtract(filtered_variants, snps)
        else if "indels" in variants:
            indels = shared_chunk[[shared_chunk["ALT"] == "NA"]]
            for i in range(-5, 5):
                indels_copy = indels.copy()
                indels_copy["POS"] = indels_copy["POS"] + i
                filtered_variants = df_subtract(
                    filtered_variants, indels_copy[["CHROM", "POS"]]
                )

    filtered_variants = filtered_variants[filtered_variants["ALT"] != "*"]
    filtered_variants.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    filter_variants()

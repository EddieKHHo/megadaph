#!/usr/bin/env python3
"""Filter variants based on depth and uniqueness.

To pass this filtering step, variant positions must:
    1. Have min_depth < depth < max_depth across all samples
    2. Show no evidence for being present across multiple samples.
    3. Be a novel heteroygous variant (not loss of heterozygosity)

Variant tables are read into memory directly so only potential variant sites
should be provided to this step.
"""
import sys

import click
from fmbiopy.df import df_subtract
from numpy import int64
from pandas import read_csv


extreme_depth = "output/find_extreme_depth_sites/GA.tsv"
shared_alleles = "output/find_shared_alleles/GA.tsv"
variants = read_csv("output/exclude_nonvariants/GA3.tsv", delim_whitespace=True)


@click.command()
@click.option(
    "--extreme-depth",
    type=click.File("r"),
    help="File containing list of sites with extreme read depth",
)
@click.option(
    "--shared_alleles",
    type=click.File("r"),
    help=(
        "File containing list of potential variants which are found in "
        "multiple samples"
    ),
)
@click.argument(
    "vcf", nargs=1, type=click.File("r"), help="VCF file to filter."
)
def filter_variants(extreme_depth, shared_alleles, variants):
    filtered_variants = read_csv(variants, sep="\t")
    for depth_chunk in read_csv(
        extreme_depth,
        delim_whitespace=True,
        chunksize=100000,
        header=None,
        names=["CHROM", "POS"],
    ):
        filtered_variants = df_subtract(filtered_variants, depth_chunk)

    for shared_chunk in read_csv(
        shared_alleles,
        delim_whitespace=True,
        chunksize=100000,
        header=None,
        names=["CHROM", "POS", "ALT"],
    ):
        filtered_variants = df_subtract(filtered_variants, shared_chunk)
    filtered_variants.to_csv(sys.stdout, sep="\t", index=False)

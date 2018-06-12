#!/usr/bin/env python3
"""Filter variants adjacent to indels in other samples."""
import sys

import click
from fmbiopy.df import df_subtract
from pandas import read_csv


def filter_noisy_regions(target_variants, base_range, other_indels):
    for indel_table in other_indels:
        for i in range(-base_range, base_range):
            indel_table_copy = indel_table.copy()
            indel_table_copy["POS"] = indel_table["POS"] + i
            target_variants = df_subtract(
                target_variants, indel_table_copy[["CHROM", "POS"]]
            )
    return target_variants


@click.command()
@click.option(
    "--target", help="Target tsv file containing variants to be filtered."
)
@click.option(
    "--nbases",
    type=int,
    help=("Number of bases surrounding variant to check in other samples."),
)
@click.argument("other_sample_indels", nargs=-1)
def cli(target, nbases, other_sample_indels):
    """Filter variants adjacent to indels in other samples."""
    other_sample_indel_tables = [
        read_csv(f, sep="\t") for f in other_sample_indels
    ]
    target_variants = read_csv(target, sep="\t")
    filtered_variants = filter_noisy_regions(
        target_variants, nbases, other_sample_indel_tables
    )
    filtered_variants.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    cli()

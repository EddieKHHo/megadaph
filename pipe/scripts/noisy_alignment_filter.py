#!/usr/bin/env python3
"""Filter variants adjacent to indels in other samples."""
import sys

import click
from pandas import read_csv


def other_samples_noisy(indel_row, base_range, other_indels):
    """Check if a SV row from target table is noisy in other samples."""
    for indel_table in other_indels:
        in_position_range = indel_table["POS"].between(
            indel_row["POS"] - base_range, indel_row["POS"] + base_range
        )
        on_correct_scaffold = indel_table["CHROM"] == indel_row["CHROM"]
        if sum(in_position_range & on_correct_scaffold):
            return True
    return False


def filter_noisy_regions(target_indels, base_range, other_indels):
    sv_is_shared = target_indels.apply(
        other_samples_noisy,
        axis=1,
        base_range=base_range,
        other_indels=other_indels,
    )
    return target_indels[~sv_is_shared]


@click.command()
@click.option(
    "--target",
    help="Target tsv file containing variants to be filtered.",
)
@click.option(
    "--nbases",
    type=int,
    help=(
        "Number of bases surrounding variant to check in other samples."
    ),
)
@click.argument("other_sample_indels", nargs=-1)
def cli(target, nbases, other_sample_indels):
    """Filter structural variants in which the alignments are noisy in other
    samples.

    This is to remove structural variants which are only called in one sample,
    but are evidenced by misalignments in other samples.
    """
    other_sample_indel_tables = [
        read_csv(f, sep="\t") for f in other_sample_indels
    ]
    target_indels_table = read_csv(target, sep="\t")
    filtered_svs = filter_noisy_regions(
        target_indels_table, nbases, other_sample_indel_tables
    )
    filtered_svs.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    cli()

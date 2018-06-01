#!/usr/bin/env python3
"""Find sites which have extreme read depth across samples.

Takes a set of readcount files produced by mpileup2readcounts and finds sites in
which the median read depth is not between `min_depth` and `max_depth`. Outputs
tsv file with 2 columns: [Chromosome/Contig, Position].
"""
import sys

import click
from pandas import DataFrame, read_csv

from fmbiopy.readcounts import destrand_counts


def find_extreme_depth(per_sample_counts, min_depth, max_depth):
    """Find positions with median coverage above `max_depth` or below
    `min_depth`.
    """
    per_sample_counts = [destrand_counts(x) for x in per_sample_counts]
    per_sample_coverage = [counts["depth"] for counts in per_sample_counts]
    medians = DataFrame(per_sample_coverage).median(axis=0)
    outside_range = (medians < min_depth) | (medians > max_depth)
    extreme_depth_sites = per_sample_counts[0].loc[outside_range][
        ["chr", "pos"]
    ]
    return extreme_depth_sites


@click.command()
@click.option("--max-depth", type=int, help="Maximum read depth")
@click.option("--min-depth", type=int, help="Minimum read depth")
@click.argument(
    "pileups", type=click.Path(exists=True, dir_okay=False), nargs=-1
)
def find_extreme_depth_sites(max_depth, min_depth, pileups):
    readcount_iter = [read_csv(f, sep="\t", chunksize=100000) for f in pileups]
    for chunks in zip(*readcount_iter):
        chunks = [chunk.rename(columns={'loc': 'pos'}) for chunk in chunks]
        extreme_depth_sites = find_extreme_depth(chunks, min_depth, max_depth)
        extreme_depth_sites.to_string(sys.stdout, header=False, index=False)
        sys.stdout.write("\n")


if __name__ == "__main__":
    find_extreme_depth_sites()

#!/usr/bin/env python
"""Filter variants for which the mean coverage across samples isn't within a
range.
"""
import click
from fmbiopy.io import copy_header
from pandas import DataFrame, read_csv
from plumbum.cmd import wc

from variant_table import VariantTable


def depth_filter(variants, min_depth, max_depth):
    """Remove variants with coverage above `max_depth` or below `min_depth`."""
    medians = DataFrame(variants.coverage).median(axis=0)
    in_range = (medians >= min_depth) & (medians <= max_depth)
    variants.df = variants.df.loc[in_range].reindex()
    return variants


def count_num_variants(tsv):
    """Count number of variants in a tsv file."""
    return int(wc("-l", tsv).rstrip()) - 1


@click.command()
@click.option("--max-depth", type=int, help="Maximum read depth")
@click.option("--min-depth", type=int, help="Minimum read depth")
@click.option(
    "--summary", type=str, help="Output text file with number of filtered sites"
)
@click.option("--output", type=str, help="Output tsv file")
@click.argument("filename", nargs=1)
def depth_filter_variants(max_depth, min_depth, summary, output, filename):
    num_filtered_variants = 0
    copy_header(filename, output)
    with open(output, "a") as outfile:
        for chunk in read_csv(filename, sep="\t", chunksize=100000):
            variants = VariantTable(chunk)
            num_variants_prefilter = len(variants)
            depth_filter(variants, min_depth, max_depth)
            num_variants_postfilter = len(variants)
            num_filtered_variants += (
                num_variants_prefilter - num_variants_postfilter
            )
            variants.df.to_csv(outfile, sep="\t", header=False, index=False)
    with open(summary, "w") as f:
        f.write(str(num_filtered_variants) + "\n")


if __name__ == "__main__":
    depth_filter_variants()

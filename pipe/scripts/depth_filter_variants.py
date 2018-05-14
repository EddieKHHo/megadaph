#!/usr/bin/env python
"""Find variants which are unique to sample from a multisample gatk variant
table.
"""
import click
from fmbiopy.io import write_table
from pandas import read_csv

from variant_table import VariantTable


def depth_filter(variants, min_depth, max_depth):
    in_range = [(depth > min_depth - 1) & (depth < max_depth + 1)
                for depth in variants.coverage]
    all_samples_in_range = sum(in_range) == variants.nsamples
    variants.df = variants.df.loc[all_samples_in_range].reindex()
    return variants


@click.command()
@click.option('--max-depth', type=int, help='Maximum read depth')
@click.option('--min-depth', type=int, help='Minimum read depth')
@click.option('--summary', type=str,
              help='Output text file with number of filtered sites')
@click.option('--output', type=str, help='Output tsv file')
@click.argument('filename', nargs=1)
def depth_filter_variants(max_depth, min_depth, summary, output, filename):
    df = read_csv(filename, sep='\t')
    variants = VariantTable(df)
    num_variants_prefilter = len(variants)
    depth_filter(variants, min_depth, max_depth)
    num_filtered_variants = num_variants_prefilter - len(variants)
    write_table(variants.df, output, sep='\t')
    with open(summary, 'w') as f:
        f.write(str(num_filtered_variants) + '\n')


if __name__ == '__main__':
    depth_filter_variants()

#!/usr/bin/env python
"""
Find variants which are unique to sample from a multisample gatk variant
table.
"""
import click
from fmbiopy.io import write_table
import pandas as pd


def count_alt_alleles(df):
    """Count number of alternative alleles per sample in gatk variant table."""
    adcols = [x for x in df if '.AD' in x]
    alt_allele_counts = [df[col].str.split(',', 1).str[0].astype(int)
                         for col in adcols]
    return alt_allele_counts


def get_unique_variants(df, max_other_counts, min_self_counts):
    altcounts = count_alt_alleles(df)
    nsamples_above_max = sum([count > max_other_counts for count in altcounts])
    nsamples_above_min = sum([count > min_self_counts for count in altcounts])
    site_is_unique = (nsamples_above_max == 1) & (nsamples_above_min == 1)
    unique_variants = df[site_is_unique].reindex()
    return unique_variants


@click.command()
@click.option('--max-other-counts', type=int,
              help='Max number of reads supporting allele in other samples.')
@click.option('--min-self-counts', type=int,
              help=('Minimum number of reads supporting allele in mutant'
                    'sample.'))
@click.option('--output', type=str, help=('Output file'))
@click.argument('filename', nargs=1)
def filter_shared_alleles(max_other_counts, min_self_counts, output, filename):
    df = pd.read_csv(filename, sep='\t')
    unique_variants = get_unique_variants(df, max_other_counts, min_self_counts)
    write_table(unique_variants, output, sep='\t')


if __name__ == '__main__':
    filter_shared_alleles()

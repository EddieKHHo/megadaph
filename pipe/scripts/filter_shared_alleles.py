#!/usr/bin/env python
"""Find variants which are unique to sample from a multisample gatk variant
table.
"""
import os

import click
from fmbiopy.io import write_table
from pandas import (
    DataFrame,
    read_csv,
)

from variant_table import VariantTable


def filter_variants(variants, cutoff):
    """Remove shared variants from a VariantTable

    Returns
    -------
    List[DataFrame]
        A DataFrame containing sites which passed filtering for each sample

    """
    counts = variants.allele_counts
    freqs = variants.allele_freq
    passing = [DataFrame() for _ in range(variants.nsamples)]
    num_alt_alleles_per_site = len(list(counts[0]))
    for allele_index in range(num_alt_alleles_per_site):
        # Loop through alleles
        nonzero_af = [freq[allele_index] > cutoff for freq in freqs]
        zero_af = [count[allele_index] == 0 for count in counts]
        for sample_index in range(variants.nsamples):
            # Loop through samples
            is_unique_hom = ((zero_af[sample_index]) &
                             (sum(nonzero_af) == variants.nsamples - 1))
            is_unique_het = ((sum(zero_af) == variants.nsamples - 1) &
                             (nonzero_af[sample_index]))
            is_unique_variant = is_unique_het | is_unique_hom
            if sum(is_unique_variant):
                passing[sample_index] = passing[sample_index].append(
                    variants.df.loc[is_unique_variant.values],
                    ignore_index=True)
    return passing


@click.command()
@click.option('--het-cutoff', type=float,
              help='Minimum allele frequency to be considered heterozygous')
@click.option('--outdir', type=str, help=('Output directory'))
@click.argument('filename', nargs=1)
def filter_shared_alleles(het_cutoff, outdir, filename):
    df = read_csv(filename, sep='\t')
    variants = VariantTable(df)
    unique_variants = filter_variants(variants, het_cutoff)
    try:
        os.mkdir(outdir)
    except OSError:
        pass
    for sample, sample_data in zip(variants.samples, unique_variants):
        outfile = os.path.join(outdir, sample + '.tsv')
        write_table(sample_data, outfile, sep='\t')


if __name__ == '__main__':
    filter_shared_alleles()

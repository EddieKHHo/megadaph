#!/usr/bin/env python3
"""Find variants which are unique to sample from mpileup2readcounts output."""
import sys

import click
from fmbiopy.readcounts import (
    destrand_counts,
    find_shared_bases,
    find_shared_indels,
    READCOUNT_DTYPE,
)
from pandas import DataFrame, read_csv

pileups = [
    "/home/fen-arch/fmacrae/megadaph.private/megadaph/pipe/output/produce_pileup/FA10.head.pileup",
    "/home/fen-arch/fmacrae/megadaph.private/megadaph/pipe/output/produce_pileup/FA7.head.pileup",
    "/home/fen-arch/fmacrae/megadaph.private/megadaph/pipe/output/produce_pileup/FAEC1.head.pileup",
]

from pandas import read_csv

readcount_iter = [
    read_csv(pileup, sep="\t", chunksize=10000) for pileup in pileups
]

per_sample_counts = next(zip(*readcount_iter))


def nrow(df):
    """Return number of rows in a DataFrame."""
    return df.shape[0]


def prepend_ref_base(indels, ref_bases):
    """Prepend reference base to indels to match vcf formatting."""
    indels_copy = indels.copy()
    for index, insertion_list in indels.iteritems():
        ref_base = ref_bases.loc[index]
        indels_copy.loc[index] = [ref_base + ins for ins in insertion_list]
    return indels_copy


def expand_list_column(column):
    """Expand a column of lists to multiple columns."""
    expanded = DataFrame(column.values.tolist())
    expanded.index = column.index
    return expanded


def find_shared_indel_pos(per_sample_counts, indel_type):
    """Find indel positions which are shared between samples.

    Parameters
    ----------
    per_sample_counts: List[DataFrame]
        List of unstranded count tables. One for each sample.
    indel_type: str
        One of ["Insertion", "Deletion"]

    Returns
    -------
    DataFrame
        DataFrame with 3 columns: ["chr", "pos", "alt"]. If multiple shared
        shared indels are found at a site, multiple rows will be returned.

    """
    shared_indel_pos = DataFrame()
    ref_bases = per_sample_counts[0]["ref"]
    # Build a dataframe of shared indels with row indices from
    # `per_sample_counts`. If multiple shared indels are present at a position,
    # split them across multiple columns.
    shared_indels = (
        find_shared_indels(per_sample_counts, indel_type)
        .pipe(prepend_ref_base, ref_bases)
        .pipe(expand_list_column)
    )
    if nrow(shared_indels) > 0:
        # Loop through indel columns, appending each to table of discovered
        # positions
        for column_index in shared_indels:
            indel_column = shared_indels[column_index].dropna()
            sequence_ids = per_sample_counts[0].loc[indel_column.index, "chr"]
            positions = per_sample_counts[0].loc[indel_column.index, "pos"]

            if indel_type == "Insertion":
                ref = ref_bases[indel_column.index]
                alt = indel_column
            elif indel_type == "Deletion":
                # Deleted bases are added to the ref column in vcf files.
                ref = indel_column
                alt = ref_bases[indel_column.index]

            shared_pos = DataFrame(
                {
                    "chr": sequence_ids.values,
                    "pos": positions,
                    "ref": ref,
                    "alt": alt,
                }
            )
            shared_indel_pos = shared_indel_pos.append(shared_pos)
    return shared_indel_pos


def find_shared_base_pos(per_sample_counts, base):
    """Find alternative bases which are shared across samples.

    Parameters
    ----------
    per_sample_counts: List[DataFrame]
        List of unstranded count tables. One for each sample.
    base: str
        One of ["A", "C", "T", "G"]

    Returns
    -------
    DataFrame
        DataFrame with 3 columns: ["chr", "pos", "alt"].

    """
    ref_bases = per_sample_counts[0]["ref"]
    is_shared_base = find_shared_bases(per_sample_counts, base) & (
        ref_bases != base
    )
    if sum(is_shared_base):
        sequence_ids = per_sample_counts[0].loc[is_shared_base]["chr"]
        positions = per_sample_counts[0].loc[is_shared_base]["pos"]
        shared_base_pos = DataFrame(
            {"chr": sequence_ids, "pos": positions, "alt": base}
        )
    else:
        shared_base_pos = DataFrame()
    return shared_base_pos


def find_shared_allele_pos(per_sample_counts):
    """Find sites which show evidence for the same allele in multiple samples.

    Parameters
    ----------
    per_sample_counts: List[DataFrame]
        List of unstranded count tables. One for each sample.

    Returns
    -------
    DataFrame

    """
    per_sample_counts = [
        destrand_counts(counts) for counts in per_sample_counts
    ]
    shared_alleles = DataFrame()
    for nuc in ["A", "C", "G", "T"]:
        shared_base_pos = find_shared_base_pos(per_sample_counts, nuc)
        shared_alleles = shared_alleles.append(
            shared_base_pos, ignore_index=True, sort=False
        )
    for indel in ["Insertion", "Deletion"]:
        shared_indel_pos = find_shared_indel_pos(per_sample_counts, indel)
        shared_alleles = shared_alleles.append(
            shared_indel_pos, ignore_index=True, sort=False
        )
    return shared_alleles


@click.command()
@click.argument("pileups", click.Path(exists=True, dir_okay=False), nargs=-1)
def find_shared_alleles(pileups):
    readcount_iter = [
        read_csv(f, sep="\t", chunksize=100000, dtype=READCOUNT_DTYPE)
        for f in pileups
    ]
    for chunks in zip(*readcount_iter):
        # Rename "loc" column so that it doesn't cause issues with loc indexing
        chunks = [chunk.rename(columns={"loc": "pos"}) for chunk in chunks]
        shared_allele_pos = find_shared_allele_pos(chunks)
        if shared_allele_pos.shape[0] > 0:
            shared_allele_pos.to_string(sys.stdout, header=False, index=False)
            sys.stdout.write("\n")


if __name__ == "__main__":
    find_shared_alleles()

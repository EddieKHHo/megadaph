#!/usr/bin/env python
"""Filter variants for which the mean coverage across samples isn't within a
range.
"""
from contextlib import ExitStack
import os

import click
from numpy import NaN
from pandas import DataFrame, isnull, read_csv
from plumbum.cmd import wc


def combine_indel_counts(indel_cell):
    """Combine indel counts from a cell in a count table."""
    if isnull(indel_cell):
        return NaN
    indel_counts = dict()
    for indel in indel_cell.split("|"):
        indel_count, indel_type = tuple(indel.split(":"))
        try:
            indel_counts[indel_type.upper()] += int(indel_count)
        except KeyError:
            indel_counts[indel_type.upper()] = int(indel_count)
    return "|".join([":".join([str(v), k]) for k, v in indel_counts.items()])


def combine_strand_counts(counts):
    """Combine all stranded counts into unstranded counts."""
    combined_counts = counts[["chr", "loc", "ref", "depth"]].copy()
    for i, nuc in enumerate(["A", "C", "T", "G"]):
        combined_counts[nuc] = counts[nuc] + counts[nuc.lower()]
    for indel_type in ["Insertion", "Deletion"]:
        combined_counts[indel_type] = counts[indel_type].apply(
            combine_indel_counts
        )
    return combined_counts


def depth_filter(per_sample_counts, min_depth, max_depth):
    """Remove variants with median coverage above `max_depth` or below
    `min_depth`.
    """
    per_sample_counts = [combine_strand_counts(x) for x in per_sample_counts]
    per_sample_coverage = [counts["depth"] for counts in per_sample_counts]
    medians = DataFrame(per_sample_coverage).median(axis=0)
    in_range = (medians >= min_depth) & (medians <= max_depth)
    filtered_counts = [
        sample_counts.loc[in_range].reindex()
        for sample_counts in per_sample_counts
    ]
    return filtered_counts


def count_num_variants(tsv):
    """Count number of variants in a tsv file."""
    return int(wc("-l", tsv).rstrip()) - 1


def write_header(output_file):
    """Write column headers to the output file."""
    output_file.write(
        "chr\tloc\tref\tdepth\tA\tC\tT\tG\tInsertion\tDeletion\t\n"
    )


@click.command()
@click.option("--max-depth", type=int, help="Maximum read depth")
@click.option("--min-depth", type=int, help="Minimum read depth")
@click.option(
    "--summary", type=str, help="Output text file with number of filtered sites"
)
@click.option("--outdir", type=str, help="Output directory")
@click.argument("pileups", nargs=-1)
def depth_filter_readcounts(max_depth, min_depth, summary, outdir, pileups):
    num_filtered_variants = 0
    output_filenames = [
        os.path.join(outdir, os.path.basename(f)) for f in pileups
    ]
    for f in output_filenames:
        if os.path.exists(f):
            os.unlink(f)

    with ExitStack() as stack:
        output_files = [
            stack.enter_context(open(fname, "a")) for fname in output_filenames
        ]
        for f in output_files:
            write_header(f)

        readcount_iter = [
            read_csv(f, sep="\t", chunksize=100000) for f in pileups
        ]
        for chunks in zip(*readcount_iter):
            num_variants_prefilter = len(chunks[0])
            all_filtered_counts = depth_filter(chunks, min_depth, max_depth)
            num_variants_postfilter = len(all_filtered_counts[0])
            num_filtered_variants += (
                num_variants_prefilter - num_variants_postfilter
            )
            for counts, output_file in zip(all_filtered_counts, output_files):
                counts.to_csv(
                    output_file,
                    sep="\t",
                    header=False,
                    index=False,
                    na_rep="NaN",
                )

    with open(summary, "w") as f:
        f.write(str(num_filtered_variants) + "\n")


if __name__ == "__main__":
    depth_filter_readcounts()

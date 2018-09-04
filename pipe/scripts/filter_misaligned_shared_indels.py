#!/usr/bin/env python3
"""Remove indels which are falsely aligned as runs of single base mismatches in
other samples."""
# This script is awful
import sys

import click

from pandas import read_csv
from plumbum.cmd import rg


BASE_IDX = {
    "A": [4, 8],
    "T": [5, 9],
    "C": [6, 10],
    "G": [7, 11]
}


def is_insertion(indel):
    return len(indel["REF"]) < len(indel["ALT"])


def get_ref_base(chrom, pos, pileup):
    exit_code, stdout, stderr = rg.run(
        ["\s".join([chrom, str(pos)]), pileup],
        retcode=None
    )
    if exit_code != 0:
        return 0
    return stdout.split()[2]


def check_indel(indel, pileups):
    if is_insertion(indel):
        bad_bases = indel["ALT"][1:]
        indel_length = len(indel["ALT"])
    else:
        indel_length = len(indel["REF"])
        bad_base_start = indel["POS"] + indel_length
        bad_base_pos = range(bad_base_start, bad_base_start + indel_length)
        bad_bases = [
            get_ref_base(indel["CHROM"], pos, pileups[0])
            for pos in bad_base_pos
        ]
        if any([base == 0 for base in bad_bases]):
            return True
    adj_base = int(indel["POS"] + 1)
    adj_base_pos = range(adj_base, adj_base + indel_length - 1)
    for bad, adj in zip(bad_bases, adj_base_pos):
        if get_ref_base(indel["CHROM"], adj, pileups[0]) != bad:
            for pileup in pileups:
                counts = rg.run(
                    ["\s".join([indel["CHROM"], str(adj)]),
                    pileup], retcode=None
                )[1].split()
                if counts:
                    target_base_counts = (
                        int(counts[BASE_IDX[bad][0]]) +
                        int(counts[BASE_IDX[bad][1]])
                    )
                    if target_base_counts > 0:
                        return False
    return True


def filter_indels(pileups, indels):
    if len(indels):
        passing = indels.apply(check_indel, axis=1, pileups=pileups)
        filtered_indels = indels[passing]
    else:
        filtered_indels = indels
    filtered_indels.to_csv(sys.stdout, sep="\t", index=False)


@click.command()
@click.option("--indels", help="TSV file containing indels.")
@click.argument("pileups", nargs=-1)
def cli(pileups, indels):
    filter_indels(pileups, read_csv(indels, sep="\t"))


if __name__ == "__main__":
    cli()

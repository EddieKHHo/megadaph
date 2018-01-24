#!/usr/bin/env python
"""Align reads with bowtie2 and sort
    print(sys.argv)

This is a wrapper around bowtie2 which handles the indexing step before
alignment, and converts the output sam to a sorted and indexed bam file.

Output is sent to stdout. -S argument should not be passed since the stdout
of Bowtie2 needs to be piped to samtools.

Requires: bowtie2, samtools

Usage:
  align_and_sort.py [args]*
  align_and_sort.py (-h | --help)

Options:
  -h --help                 Show this screen
"""
from __future__ import print_function
import sys

from plumbum import (
        BG,
        FG,
        local,
        )
from plumbum.cmd import (
        bowtie2,
        samtools,
        )

def _get_reference(args):
    """Get the path to the reference seq from a list of bowtie2 arguments"""
    x_index = args.index('-x')
    ref_index = local.path(args[x_index + 1])

    reference = local.path(ref_index + '.fa')
    if not reference.exists():
        reference = local.path(ref_index + '.fasta')
        if not reference.exists():
            raise ValueError('Fasta not found')
    return reference


def align_and_sort(args):
    """Main function"""
    ref = _get_reference(args)
    fai_indices = local.path(ref + '.fai')
    if not fai_indices:
        print('Indexing fasta with samtools', sys.stderr)
        samtools_index = samtools('faidx', ref) & BG
        samtools_index.wait()

    bowtie2_index = ref.with_suffix('.1.bt2')
    if not bowtie2_index.exists():
        print('Indexing fasta with bowtie2')
        prefix = ref.with_suffix('')
        build = local['bowtie2-build']
        p = build[ref, prefix] & BG
        p.wait()

    align_pipe = (
            bowtie2.__getitem__(args) |
            samtools['view', '-bhS'] |
            samtools['sort'])
    p = align_pipe & BG(stdout = sys.stdout, stderr = sys.stderr)
    p.wait()

if __name__ == '__main__':
    if sys.argv[1] in ['-h', '--help']:
        print(__doc__)
        sys.exit(0)
    align_and_sort(sys.argv[1:])

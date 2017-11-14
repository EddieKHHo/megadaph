#!/usr/bin/env python
"""Rename the original Megadaph files, using a simpler naming convention

Example:

lane8-s046-indexD711-D501-TCTCGCGC-AGGCTATA-A10_73_FA7_S46_L008_R2_001.fq.gz

-> FA7.R2.fq.gz
"""
from sys import argv
from typing import Sequence

from fmbiopy.fmpaths import as_paths

def rename_megadaph_files(files: Sequence[str]):
    """Rename the original Megadaph files, using a simpler naming convention"""
    # Convert to a list of paths
    paths = as_paths(files)
    for path in paths:
        suffix = ''.join(path.suffixes)
        name = path.name
        unscore_split = name.split('_')
        sample = unscore_split[2]
        pair = unscore_split[5]
        new_name = ''.join([sample, '.', pair, suffix])
        path.rename(new_name)


if __name__ == "__main__":
    if len(argv) > 2:
        inputs = argv[1:]
    else:
        inputs = [argv[1]]
    rename_megadaph_files(inputs)

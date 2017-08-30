#!/usr/bin/env python

""" Main pipeline script for Megadaph assembly decontamination 

Assumes sort order 

Usage:
    megadaph_decontam_pipeline.py -r <read_dir> <assembly>...

Options:
    -h --help   Show this screen
    -r          Directory where reads are located

""" 

import os, sys
import fen_util
from docopt import docopt
from modules.alignment import alignment
from modules.parse_arguments import parse_arguments

if __name__ == '__main__':
    args = docopt(__doc__)
    fwd_reads, rev_reads, assemblies = parse_arguments(args)
    arg_files = get_arg_files(args)

"""Produce a series of blobplots across varius taxonomic ranks

Requires 5 threads

Usage:
  produce_blobplots.py -o DIR -d BLOBTOOLS_DATABASE

Options:
  -o, --outdir=DIR              Output directory
  -d, --db=BLOBTOOLS_DATABASE   Blobtools database json
"""
from docopt import docopt
from plumbum import local
from plumbum.cmd import blobtools


RANKS = ["family", "genus", "family", "order", "phylum", "species"]


def main(db, outdir):
    ps = []
    for rank in RANKS:
        rank_outdir = local.path(outdir) / rank
        ps.append(
            blobtools.popen('blobplot', '-p', '20', '-r', rank, '-i', db, '-o',
                            rank_outdir))
    for p in ps:
        p.communicate()


if __name__ == "__main__":
    opts = docopt(__doc__)
    main(opts['--db'], opts['--outdir'])

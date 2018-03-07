#!/usr/bin/env python
'''Remove mitochondrial contigs from the megadaph assemblies'''
import csv
import logging
from logging import info as log
import os
import sys

from Bio import SeqIO
from plumbum import local
import click

logging.basicConfig(level=os.environ.get('LOGLEVEL', 'INFO'))


@click.command()
@click.option('--blast', help='Blast tsv file with outfmt=6',
              type=click.Path(exists=True, dir_okay=False))
@click.option('--assembly', help='Assembly file (fasta)',
              type=click.Path(exists=True, dir_okay=False))
@click.option('--outfile', help='Output filename', type=click.Path())
@click.option('--min_blast_length', default=2000,
              help=('Minimum BLAST alignment length for a scaffold to be '
                    'considered mitochondrial'))
@click.option('--max_scaffold_length', default=20000,
              help=('Maximum scaffold size for a scaffold to be considered '
                    'mitochondrial'))
def main(blast, assembly, outfile, min_blast_length, max_scaffold_length):
    '''Remove mitochondrial contigs from the megadaph assemblies'''
    blast = local.path(blast)
    assembly = local.path(assembly)
    outfile = local.path(outfile)
    log('Finding candidates mitochondrial sequences...')
    mt_candidates = get_mt_seqids(blast, min_blast_length)

    if mt_candidates:
        log('Found the following candidates: %s', mt_candidates)
        log('Removing candidates from assembly...')
        remove_mt_from_assembly(mt_candidates, assembly, max_scaffold_length,
                                outfile)
    else:
        log('No candidates found, copying input to output...')
        assembly.copy(outfile)
    log('Shutting down...')
    sys.exit()


def remove_mt_from_assembly(mt_candidates, assembly, max_scaffold_length,
                            outfile):
    '''Remove mitochondrial scaffolds from the assembly'''
    with assembly.open('r') as inp, outfile.open('w') as out:
        for record in SeqIO.parse(inp, 'fasta'):
            if record.id in mt_candidates:
                sequence_length = len(record)
                if sequence_length < max_scaffold_length:
                    log('Removed %s from assembly (length = %s)', record.id,
                        sequence_length)
                else:
                    log('Not removing %s from assembly because its too long '
                        '(length = %s)', record.id, sequence_length)
                    SeqIO.write(record, out, 'fasta')
            else:
                SeqIO.write(record, out, 'fasta')


def get_mt_seqids(blast, min_blast_length):
    '''Get candidate mitochondrial sequence IDs from the BLAST file'''
    seq_ids = []
    with blast.open('r') as f:
        tsv = csv.reader(f, delimiter='\t')
        for row in tsv:
            alignment_length = int(row[3])
            seq_id = row[1]
            if alignment_length >= min_blast_length:
                seq_ids.append(seq_id)
    return seq_ids

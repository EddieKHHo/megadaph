#!/usr/bin/env/ python
import os
import fmbiopy.fmbiopy.fen_util as fen_util
from fmbiopy.fmbiopy.run_bowtie2 import run_bowtie2
from fmbiopy.fmbiopy.run_index_fasta import run_index_fasta

## Align reads to their corresponding assembly using bowtie2, and convert to 
##.bam
## Param
##   fwd_reads, rev_reads   Fastq sequences, gzip supported
def alignment(fwd_reads, rev_reads, references, out_bam_paths, threads=None, 
              maxins=None, no_discordant=None, no_mixed=None, preset=None):
    
    # Check arguments
    fen_util.check_all_exist(fwd_reads + rev_reads + references)
    fen_util.check_all_suffixes(fwd_read + rev_reads, [".fq", ".fq.gz",
        ".fastq", ".fastq.gz"])

    indices = run_index_fasta(references)

    out_sam_paths = []
    for p in out_bam_paths:
        sam_path = fen_util.replace_suffix(p, "bam", "sam")
        out_sam_paths.append(sam_path)

    sam_files = run_bowtie2(fwd_reads, rev_reads, indices, out_sam_paths,
                            threads, maxins, no_discordant, no_mixed, preset)

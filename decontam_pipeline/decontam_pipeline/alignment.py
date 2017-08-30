#!/usr/bin/env/ python
import fen_util
import os
from run_bowtie2 import run_bowtie2
from run_index_fasta import run_index_fasta

def alignment(fwd_reads, rev_reads, references, out_bam_paths, threads=None, 
              maxins=None, no_discordant=None, no_mixed=None, preset=None):
    
    ref_dir = os.path.dirname(references[0])
    with fen_util.working_directory(ref_dir):
        indices = run_index_fasta(references)

    out_sam_paths = []
    for p in out_bam_paths:
        sam_path = fen_util.replace_suffix(p, "bam", "sam")
        out_sam_paths.append(sam_path)

    sam_files = run_bowtie2(fwd_reads, rev_reads, indices, out_sam_paths,
                            threads, maxins, no_discordant, no_mixed, preset)

#!/usr/bin/env python
"""Snakemake script for running Haplomerger2"""
from shutil import move
from subprocess import CalledProcessError
import sys

from plumbum import (
    FG,
    local,
)
from plumbum.cmd import (
    gzip,
    windowmasker,
)
from snakemake import shell

assembly = local.path(snakemake.input[0])
log = local.path(snakemake.log[0])
try:
    threads = snakemake.params.threads
except:
    threads = "1"


# Copy the template project directory
template.copy(projdir)

# Create a masking library
masklib = projdir / "masklib.ustat"
windowmasker["-checkdup", "true", "-mk_counts", "-in", assembly, "-out",
             masklib, "-mem", "10000"] & FG

# Soft mask assembly
masked_assembly = projdir / (snakemake.wildcards.sample + ".fa")
windowmasker["-ustat", masklib, "-in", assembly, "-out", masked_assembly,
             "-outfmt", "fasta", "-dust", "true"] & FG
gzip[masked_assembly] & FG
masked_assembly = local.path(masked_assembly + ".gz")

with local.cwd(projdir):
    try:
        shell(' '.join(["./run_all.batch", snakemake.wildcards.sample, threads,
                        "2>", str(log)]))
    except CalledProcessError:
        pass

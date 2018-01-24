#!/usr/bin/env python
"""Snakescript for softmasking the initial assemblies"""
import os

from snakemake import shell

masklib = os.path.join(snakemake.input.projdir, "masklib.ustat")
shell(' '.join(["windowmasker -checkdup true -mk_counts -in",
                snakemake.input.assembly, "-out", masklib, "-mem 10000"]))
shell(' '.join(["windowmasker -ustat", masklib, "-in", snakemake.input.assembly,
                "-outfmt fasta -dust true | gzip - >", snakemake.output[0]]))

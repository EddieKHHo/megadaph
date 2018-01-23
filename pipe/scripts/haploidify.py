#!/usr/bin/env python
"""Snakescript for remove misjoins in the initial megadaph assemblies"""
from plumbum import local
from snakemake import shell

try:
    threads = str(snakemake.threads)
except NameError:
    threads = 1

inp = local.path(snakemake.input[0])
hm2_dir = inp.dirname
sample_id = snakemake.wildcards.sample + "_A"

with local.cwd(hm2_dir):
    shell(' '.join(["./hm.batchB1.initiation_and_all_lastz", sample_id,
                    threads]))
    shell(' '.join(["./hm.batchB2.chainNet_and_netToMaf", sample_id, threads]))
    shell(' '.join(["./hm.batchB3.haplomerger", sample_id, threads]))
    shell(' '.join(["./hm.batchB4.refine_unpaired_sequences", sample_id,
                    threads]))
    shell(' '.join(["./hm.batchB5.merge_paired_and_unpaired_sequences",
                    sample_id, threads]))

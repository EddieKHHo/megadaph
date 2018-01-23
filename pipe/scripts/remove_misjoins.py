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
sample_id = snakemake.wildcards.sample

with local.cwd(hm2_dir):
    shell(' '.join(["./hm.batchA1.initiation_and_all_lastz", sample_id,
                    threads]))
    shell(' '.join(["./hm.batchA2.chainNet_and_netToMaf", sample_id, threads]))
    shell(' '.join(["./hm.batchA3.misjoin_processing", sample_id, threads]))

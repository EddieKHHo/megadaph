#!/usr/bin/env python
"""Convert bam files to blobtools COV files"""
# Snakescript
from snakemake import shell

bam_args = ' -b '.join(snakemake.input.bams)

shell("blobtools map2cov -i {snakemake.input.assembly} "
      "--output {snakemake.params.output_prefix} " +
      bam_args + " > {snakemake.log} 2>&1")

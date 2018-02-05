#!/usr/bin/env python
"""Initialize a blobtools database and move the resultant coverage files"""
# Snakescript
from plumbum import local
from snakemake import shell

bam_args = ' -b '.join(snakemake.input.bams)
output_prefix = (local.path("output/init_blobtools_database1") /
                 snakemake.wildcards.genotype)

shell("blobtools create -i {snakemake.input.assembly} --title "
      "{snakemake.wildcards.genotype} --out " + output_prefix + " --hitsfile "
      "{snakemake.input.blast_hits} --hitsfile {snakemake.input.diamond_hits} "
      "--nodes {snakemake.config[nodes_dmp]} --names "
      "{snakemake.config[names_dmp]} -b " + bam_args + " > {snakemake.log} "
      "2>&1")

# The coverage files are produced in the same directory as the database by
# default. Move them to their own directory
outdir = output_prefix.dirname
covdir = local.path(snakemake.params.covdir)
covdir.mkdir()
for x in outdir.glob("*.cov"):
    sample = x.name.split(".")[1]
    destination = covdir / (sample + ".cov")
    x.move(destination)

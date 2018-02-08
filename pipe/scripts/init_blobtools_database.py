:!/usr/bin/env python
"""Initialize a blobtools database and move the resultant coverage files"""
# Snakescript
from plumbum import local
from snakemake import shell

cov_args = ' -c '.join(snakemake.input.covs)

shell("blobtools create -i {snakemake.input.assembly} --title "
      "{snakemake.wildcards.genotype} --out " + output_prefix + " --hitsfile "
      "{snakemake.input.blast_hits} --hitsfile {snakemake.input.diamond_hits} "
      "--nodes {snakemake.config[nodes_dmp]} --names "
      "{snakemake.config[names_dmp]} -c " + cov_args + " > {snakemake.log} "
      "2>&1")

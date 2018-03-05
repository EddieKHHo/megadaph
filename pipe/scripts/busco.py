#!/usr/bin/env python
'''Snakescript for running busco on assemblies'''
from plumbum import (
    FG,
    local,
)
from plumbum.cmd import busco


outdir = local.path(snakemake.params.outdir)
outdir.delete()
outdir.mkdir()

outdir_root = outdir.dirname

local.env['BUSCO_CONFIG_FILE'] = str(local.path(snakemake.config['busco_config']))
local.env['AUGUSTUS_CONFIG_PATH'] = str(local.path(
    snakemake.config['augustus_config']))

assembly = local.path(snakemake.input)
lineage_file = local.path(snakemake.config['busco_database'])
tmp_outdir = outdir_root / ('run_' + snakemake.wildcards.genotype)
tmp_outdir.mkdir()

with local.cwd(tmp_outdir):
    busco['-i', assembly, '-o', snakemake.wildcards.genotype, "-l",
          lineage_file, '-m', 'geno', '-f'] & FG

for f in tmp_outdir.list():
    f.move(outdir)

tmp_outdir.delete()

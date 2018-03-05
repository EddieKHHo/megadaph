'''Snakescript for running busco on assemblies'''
from snakemake import shell
from plumbum import local

local.env['BUSCO_CONFIG_FILE'] = local.path(snakemake.config['busco_config'])
local.env['AUGUSTUS_CONFIG_PATH'] = local.path(
    snakemake.config['augustus_config'])

output_dir = (local.path('output') /
              ('busco' + snakemake.wildcards.iter) /
              snakemake.wildcards.genotype)
assembly = local.path(input)
lineage_file = local.path(snakemake.config['busco_database'])
busco = local['busco']
tmp_outdir = output_dir / ("run_" + snakemake.wildcards.genotype)

with local.cwd(output_dir):
    busco['-i', assembly, '-o', snakemake.wildcards.genotype, "-l",
          lineage_file, '-m', 'geno', '-f'] & FG

    for f in tmp_outdir.list():
        f.move(output_dir)

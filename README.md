Estimating the nuclear mutation rate in *Daphnia magna*
=======================================================

This repo contains the analysis pipeline for a mutation accumulation experiment *D. magna* (2017--). The pipeline is implemented in [Snakemake] ([docs]). See [pipe/Snakefile] for the main script. For a summary of the pipeline implementation see [pipe/README.md].

Dependencies:

-   https://github.com/fennerm/fmbiopy
-   https://github.com/fennerm/fen.R.util
-   snakemake
-   samtools
-   SPAdes
-   diamond blastx
-   sambamba
-   blobtools
-   picard
-   bowtie2
-   BBMap
-   blastn
-   QUAST
-   BUSCO
-   fastqc
-   multiqc
-   bioawk
-   GATK

  [Snakemake]: https://bitbucket.org/snakemake/snakemake/wiki/Home
  [docs]: https://snakemake.readthedocs.io/en/stable
  [pipe/README.md]: https://github.com/fennerm/megadaph/blob/master/pipe/README.md
  [pipe/Snakefile]: https://github.com/fennerm/megadaph/blob/master/pipe/Snakefile

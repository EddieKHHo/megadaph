#!/bin/bash
module load bedtools
module load samtools
module load perl

cd $1

###get coverage
samtools faidx mapping/$1.scaffolds.fasta
awk '{print $1, $2}' mapping/$1.scaffolds.fasta.fai >genome.txt
bedtools genomecov -ibam mapping/merged.bam -g genome.txt >genome_cov.txt
../CalcMedianCov.pl genome_cov.txt medianCov.txt

###get gc and length
../get_gc_content.pl mapping/$1.scaffolds.fasta

###combine gc, length, coverage
../combineCovGCv1.pl medianCov.txt gc_out.txt $1.gc.len.cov.txt
rm -rf genome.txt
rm -rf genome_cov.txt
rm -rf medianCov.txt
rm -rf gc_out.txt

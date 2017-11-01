#!/bin/bash

cd $1

module load bwa
module load samtools
bwa index -a is ../../$1/scaffolds.fa

bwa mem -t 8 ../../$1/scaffolds.fa ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.R1.fastq.gz ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.R2.fastq.gz >aln.pe.sam
bwa mem -t 8 ../../$1/scaffolds.fa ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.merged.fastq.gz >aln.se.sam

samtools view -bhSo aln.se.bam aln.se.sam
samtools sort -@ 8 aln.se.bam aln.se.sort

samtools view -bhSo aln.pe.bam aln.pe.sam
samtools sort -@ 8 aln.pe.bam aln.pe.sort

samtools merge -@ 8 merged.bam aln.pe.sort.bam aln.se.sort.bam

#!/bin/bash
module load samtools
module load python
 
cd $1
#awk '{if($2>30 && $2<55 && $4>5){print $0}}' $1.gc.len.cov.txt | awk '{print $1}' >list
#xargs samtools faidx mapping/$1.scaffolds.fasta < list > $1.30_55.5x.fasta
 
~/src/redundans/redundans.py -v -f $1.30_55.5x.fasta -i ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.R1.fastq.gz ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.R2.fastq.gz -l ~/megadaph/reads/cleaned_reads/Aug2017/$1/$1.merged.fastq.gz -r ../../genomes/daphnia/dmagna-v2.4-20100422-assembly.fna -o redundans -t 8


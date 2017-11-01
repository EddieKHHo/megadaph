#!/bin/bash

mkdir $1
cd $1

R1_NAME=../../../../megadaph/reads/cleaned_reads/Aug2017/$1/$1.R1.fastq.gz
R2_NAME=../../../../megadaph/reads/cleaned_reads/Aug2017/$1/$1.R2.fastq.gz
MERGE_NAME=../../../../megadaph/reads/cleaned_reads/Aug2017/$1/$1.merged.fastq.gz

../../../../src/SPAdes-3.9.1-Linux/bin/spades.py --pe1-1 $R1_NAME --pe1-2 $R2_NAME --s1 $MERGE_NAME -o $1.spades

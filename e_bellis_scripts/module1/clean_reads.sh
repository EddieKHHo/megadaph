#!/bin/bash
##need to have file of adaptor sequences in same directory

SAMPLE_ID=$1

###########PART1: copy files into new directory and check number of raw reads
mkdir $SAMPLE_ID
cd $SAMPLE_ID
cp ~/megadaph/reads/illumina_150bp/*$SAMPLE_ID* .
gunzip *gz
READ_ID=$(head -n 1 lane*R1*fastq | awk '{print substr($0,0,7)}')
grep -c "$READ_ID" *fastq

###########PART 2: trim adaptors
R1_NAME=$(ls lane*R1*fastq)
R2_NAME=$(ls lane*R2*fastq)
~/src/bbmap/bbduk.sh in=$R1_NAME in2=$R2_NAME out=adaptor_trimmed_R1.fastq out2=adaptor_trimmed_R2.fastq ref=../adaptors_illumina.fa k=23 ktrim=r mink=4 hdist=1 tpe tbo

###########PART 3: merge overlapping reads
~/src/bbmap/bbmerge.sh in1=adaptor_trimmed_R1.fastq in2=adaptor_trimmed_R2.fastq outm=merged.fastq outu1=unmerged_R1.fastq outu2=unmerged_R2.fastq t=8 ihist=hist.txt vstrict=t

###########PART 4: quality filtering
~/src/bbmap/bbduk.sh in=unmerged_R1.fastq in2=unmerged_R2.fastq out=clean.R1.fastq out2=clean.R2.fastq qtrim=rl trimq=20 minlen=50
~/src/bbmap/bbduk.sh in=merged.fastq out=clean.merged.fastq qtrim=rl trimq=20 minlen=50

###########PART 5: clean up
grep -c "$READ_ID" clean*fastq
mv clean.R1.fastq $SAMPLE_ID.R1.fastq
mv clean.R2.fastq $SAMPLE_ID.R2.fastq
mv clean.merged.fastq $SAMPLE_ID.merged.fastq
gzip $SAMPLE_ID.*.fastq

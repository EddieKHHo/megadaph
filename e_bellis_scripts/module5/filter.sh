#!/bin/bash

cd $1

cat $1.blob.txt |awk '{if($3 >= 5000 && $4 >= 33 && $4 <= 76){print $1}}'>list
xargs samtools faidx ../../$1/redundans/scaffolds.reduced.fa < list > $1.filtered.fa

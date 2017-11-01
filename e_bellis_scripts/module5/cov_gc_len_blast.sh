#!/bin/bash
module load bedtools
module load samtools
module load perl
module load bwa
module load ncbi

cd $1

###get coverage
samtools faidx ../../$1/redundans/scaffolds.reduced.fa
awk '{print $1, $2}' ../../$1/redundans/scaffolds.reduced.fa.fai >genome.txt
bedtools genomecov -ibam merged.bam -g genome.txt >genome_cov.txt
awk '{if($2>0){print $0}}' genome_cov.txt > genome_cov2.txt
../../CalcMedianCov.pl genome_cov2.txt medianCov.txt

###get gc and length
../../get_gc_content.pl ../../$1/redundans/scaffolds.reduced.fa
grep -v 'ID' gc_out.txt | awk '{GC2 = ($4+$5)/($4+$5+$6+$7)*100; print $0, GC2}' | sed -r 's/\s/\t/g' >gc_out2.txt

###combine gc, length, coverage
../../combineCovGCv2.pl medianCov.txt gc_out2.txt $1.gc.len.cov.txt
grep -v 'Contig' $1.gc.len.cov.txt | sort -k1,1 >glc.sort.txt

rm -rf genome.txt
rm -rf genome_cov.txt
rm -rf medianCov.txt
rm -rf gc_out.txt
rm -rf gc_out2.txt
rm -rf genome_cov2.txt

### blast to D. magna reference genome
blastn -task megablast -subject ../../../genomes/daphnia/dmagna-v2.4-20100422-assembly.fna -query ../../$1/redundans/scaffolds.reduced.fa -evalue 1e-6 -max_target_seqs 1 -max_hsps 1 -outfmt '6 qseqid sseqid length evalue pident qcovhsp' -out megablast_daph.txt
sort -k1,1 megablast_daph.txt >blast.sort.txt
join -a 1 -1 1 -2 1 glc.sort.txt blast.sort.txt >$1.blob.txt

rm -rf megablast_daph.txt
rm -rf blast.sort.txt
rm -rf glc.sort.txt
rm -rf $1.gc.len.cov.txt

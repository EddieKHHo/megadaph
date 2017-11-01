#!/usr/bin/perl

#Objective: this script will combine two files listing coverage and gc content

$covfile = $ARGV[0];		# file of modal coverages, one row per contig
$gcfile = $ARGV[1];		#file of gc contents
$outfile = $ARGV[2];		# name for output file, .txt format

# loop through fastq file and truncate sequences and quality scores
open (IN, "<", $covfile);
open(IN2, "<", $gcfile);
open (OUT, ">", $outfile);

print OUT "Contig\tGC\tLength\tModal_coverage\n";

%contigCov = ();

while ($line = <IN>) {
	chomp($line);
	@lineParts = split(/\t/, $line);
	$contigCov{$lineParts[0]} = $lineParts[1];
}

while ($line = <IN2>) {
	chomp($line);
	@lineParts = split(/\t/, $line);
	if ( $lineParts[0] =~ />(NODE_[0-9]+_length_[0-9]+_cov_[0-9]+.[0-9]+)/) {
		print OUT "$1\t$lineParts[1]\t$lineParts[2]\t";
		print OUT $contigCov{$1}."\n";
	}	
}

close(IN);
close(IN2);
close(OUT);

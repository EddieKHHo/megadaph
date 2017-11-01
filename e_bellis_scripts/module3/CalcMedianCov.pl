#!/usr/bin/perl

##Objective: This script will calculate median coverage of each sequence from a coverage file output by bedtools genomecov'\n"

$covfile = $ARGV[0];		# coverage file with 5 columns, .txt format
$outfile = $ARGV[1];		# name for output file, .txt format

# loop through fastq file and truncate sequences and quality scores
open (IN, "<", $covfile);
open (OUT, ">", $outfile);

%contigCovs = ();

$line = <IN>;
chomp($line);
@lineParts = split(/\t/, $line);
$prevContig = $lineParts[0];

%contigCov = ();
$contigCov{$lineParts[1]} = $lineParts[2];

while ($line = <IN>) {
	chomp($line);
	@lineParts = split(/\t/, $line);

	if($lineParts[0] =~ /$prevContig/){
		$contigCov{$lineParts[1]} = $lineParts[2];

	} else {
		@coverages = sort { $contigCov{$a} <=> $contigCov{$b} } keys %contigCov;
		$highest = $coverages[-1];
		print OUT "$prevContig\t$highest\n";
		$prevContig = $lineParts[0];

		%contigCov = ();
		$contigCov{$lineParts[1]} = $lineParts[2];
		next;
	}
}

@coverages = sort { $contigCov{$a} <=> $contigCov{$b} } keys %contigCov;
$highest = $coverages[-1];
print OUT "$prevContig\t$highest\n";

close(IN);
close(OUT);

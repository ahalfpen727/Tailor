#!/usr/bin/perl -w
use strict;

my ($knownToGeneFile, $geneExpFile, $outFile) = ($ARGV[0], $ARGV[1], $ARGV[2]);
$outFile //= "geneExpToRef.txt";

my %knownToGene;

open K2E, "<$knownToGeneFile" or die "Cannot open input file $knownToGeneFile: $!\n";
<K2E>;
while(<K2E>){
	my ($ucId, $gene) = split(/\s+/, $_);
	$knownToGene{$ucId} = $gene;
}
close K2E;

open IN, "<$geneExpFile" or die "Cannot open input file $geneExpFile: $!\n";
open OUT, ">$outFile" or die "Cannot create output file $outFile: $!\n";
my $header = <IN>;
$header = "chr\tstart\tstop\t$header";
print OUT "$header";
$, = "\t";
while(<IN>){
	#print "$_";
	my @vals = split(/\s+/, $_);
	if(defined $knownToGene{$vals[0]}){
		unshift(@vals, $knownToGene{$vals[0]});
	} else {
		unshift(@vals, "-");
	}
	my @coords = split(/:|-/, $vals[3]);
	unshift(@vals, @coords);
	#print @coords, "---$vals[3]\n";
	print OUT join("\t", @vals)."\n";
}
close IN;
close OUT;

#!/usr/bin/perl
#############################################
#This program extracts the desired sequence form a fasta file using the coordinates. 
#This program takes the geneome file in fasta format and coordinate file in tab format with the following columns: Strand, Start and End. 
# 
#Author: Yogendra Bhaskar
#
#License: GPL 3.0
#
#############################################

use strict;
use warnings;

if(!@ARGV){
	print "No FilePath Specified!\n";
	print "Usage: perl extract_seq.pl -h\n";	
	exit;
}

if($ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
	help();
	exit;
}

(@ARGV == 2) || die "Usage: perl extract_seq.pl <genome_fasta_file> <coordinate_file>\n";

sub help { print "
Usage: perl extract_seq.pl <genome_fasta_file> <coordinate_file>\n
This program takes two files as input:
1. geneome file in fasta format. 
2. coordinate file in tab format with the following columns: Strand, Start and End (readme).
Format example for coordinate file is as follows:

Strand	| Start	| End |
0     	| 150	| 200 |
1     	| 210	| 250 |
0     	| 270	| 340 |
\n";}

my $fas_file = $ARGV[0];
my $coords_file = $ARGV[1];


open INFILE1, $fas_file or die "Could not open $fas_file: $!";
my $fasta;

<INFILE1>; 

while(<INFILE1>) {
    chomp;
    $fasta .= $_;
}
close INFILE1 or die "Could not close '$fas_file'. $!";

open INFILE, $coords_file or die "Could not open $coords_file: $!";
my @read1=<INFILE>;
my $count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

my (@strand, @start, @end);
for(my $i=0; $i<$count_read1; $i++){
	my @values1= split("\t", $read1[$i]);
	$strand[$i]= $values1[0];
	$start[$i]= $values1[1];
	$end[$i]= $values1[2];
}


for(my $i=0; $i<$count_read1; $i++){
	if($strand[$i]==1){
		my $offset = ($end[$i]- ($start[$i] -10)) +1;
		my $sub= substr($fasta, ($start[$i]-10) -1, $offset);
		my $rev_comp= `rev<<<$sub| tr ATGCatgc TACGtacg|tr -d "\n"`;
		my $ext10= substr($rev_comp, (length($rev_comp) -10), 10);
		print "$strand[$i]\t$start[$i]\t$end[$i]\t$rev_comp\t$ext10\n";
	}
	elsif($strand[$i]==0){
		my $offset = ($end[$i]+10) - ($start[$i] - 1);
		my $sub = substr ($fasta, ($start[$i]-1), $offset);
		my $ext10= substr($sub, (length($sub) -10), 10);
		print "$strand[$i]\t$start[$i]\t$end[$i]\t$sub\t$ext10\n";
	}
}


close INFILE or die "Could not close '$coords_file'. $!";

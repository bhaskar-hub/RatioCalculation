#!/usr/bin/perl
#############################################
#This program calculates the maximum U content and Poly-U tails in a given sequence. 
# 
#Author: Yogendra Bhaskar
#
#License: GPL 3.0
#
#############################################

use strict;
use warnings;
use List::Util qw(min max);
# provide the extended sequences (extended_right_stem_seq)

if(!@ARGV){
	print "No FilePath Specified!\n";
	print "Usage: perl max_T_cont.pl <sequence_file>\n";
	exit;
}

if($ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
	help();
	exit;
}

sub help { print "
Usage: perl max_T_cont.pl <sequence_file>\n
This program takes a plain sequence file with each line having a single sequence (readme).
\n";}

my $sequence= $ARGV[0];
open fh1,$sequence;
#print "file opened for reading \n";
my @read1=<fh1>;
my $count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

my @seq; my $length; my $sub; my @t;
my ($max, $t3, $t4, $t5);
for(my $i=0; $i<$count_read1; $i++){
	my @values1= split("\t", $read1[$i]);
	$seq[$i]=$values1[0];
	$length=length$seq[$i];
	#print "$length\n";
	for(my $j=0; $j<$length; $j++){
		$sub= substr ($seq[$i],$j,10 );
		#calculating T content
		$t[$j]=($sub=~tr/t//);
		#print "$sub\t$t\n";
	}
	$max=max(@t);
	$t3=($seq[$i]=~s/ttt/ttt/g);
	$t4=($seq[$i]=~s/tttt/tttt/g);
	$t5=($seq[$i]=~s/ttttt/ttttt/g);
	print "$seq[$i]\t$max\t$t3\t$t4\t$t5\n";
	
}


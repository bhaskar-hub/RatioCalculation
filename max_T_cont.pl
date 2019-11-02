#!/usr/bin/perl
#use strict;
#use warnings;
use List::Util qw(min max);
# provide the extended sequences (extended_right_stem_seq)

$sequence= $ARGV[0];
open fh1,$sequence;
#print "file opened for reading \n";
@read1=<fh1>;
$count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

for($i=0; $i<$count_read1; $i++){
	@values1= split("\t", $read1[$i]);
	$seq[$i]=@values1[0];
	$length=length$seq[$i];
	#print "$length\n";
	for($j=0; $j<$length; $j++){
		$sub= substr ($seq[$i],$j,10 );
		#calculating T content
		@t[$j]=($sub=~tr/t//);
		#print "$sub\t$t\n";
	}
	$max=max(@t);
	$t3=($seq[$i]=~s/ttt/ttt/g);
	$t4=($seq[$i]=~s/tttt/tttt/g);
	$t5=($seq[$i]=~s/ttttt/ttttt/g);
	print "$seq[$i]\t$max\t$t3\t$t4\t$t5\n";
	
}


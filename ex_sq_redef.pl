#!/usr/bin/perl
#use strict;
#use warnings;

my $fas_file = $ARGV[0];
my $coords_file = $ARGV[1];


open INFILE1, $fas_file or die "Could not open $fas_file: $!";
my $fasta;

<INFILE1>; 

while(<INFILE1>) {
    chomp;
    $fasta .= $_;
}
close INFILE1 or die "Could not close '$chrom'. $!";

open INFILE, $coords_file or die "Could not open $coords_file: $!";
@read1=<INFILE>;
$count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

for($i=0; $i<$count_read1; $i++){
	@values1= split("\t", $read1[$i]);
	$strand[$i]= $values1[0];
	$start[$i]= $values1[1];
	$end[$i]= $values1[2];
}


for($i=0; $i<$count_read1; $i++){
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

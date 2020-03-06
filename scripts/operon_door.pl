#!/usr/bin/perl

use strict;
use warnings;
if(!@ARGV){
	print "No FilePath Specified!\n";
	print "Usage: perl operon_door.pl -h\n";
	exit;
}

if($ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
	help();
	exit;
}

sub help { print "This program is Perl and UNIX utilities-based tool, which parses operon annotation from DOOR database into a simple five-column (#operon, genes, strand, operon_start and operon_end) tab delimeted file. 
Usage: perl operon_door.pl <operon_file_from_door_database> 
Output file name will be \"operon_map_door\" \n";}


my $op = $ARGV[0];

system_commands($op);

sub system_commands{

	system('sed \'1d\' ' .$op.' |cut -f1,3,4,5,6| sort -k3 -g|awk \'{FS="\t"}{print $1"\t"$2"\t"$5"\t"$3"\t"$4}\'|sed \'1d\'  > operon_struc_door');
	}

	
open(fh1,"operon_struc_door");
#print "file opened for reading \n";
my @read1=<fh1>;
my $count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

my (@OP, @start, @end, @gene, @strand);
for(my $i=0; $i<$count_read1; $i++){
	my @values1= split("\t", $read1[$i]);
	$OP[$i]= $values1[0];
	$start[$i]= $values1[3];
	$end[$i]= $values1[4];
	$gene[$i]= $values1[1];
	$strand[$i]= $values1[2];
}

open(de_dup, '>de_dup_operon_struc_door');

my $b=0; my $f_sum=0; my $r_sum=0;
for (my $i=0; $i < $count_read1; $i++){
	if($OP[$i]==$OP[$i-1]){
		$gene[$i]= "$gene[$i-1],$gene[$i]";
		$start[$i]= $start[$i-1];
		$end[$i]= $end[$i];
		print de_dup "$OP[$i]\t$strand[$i]\t$gene[$i]\t$start[$i]\t$end[$i]\n";
	}
	else{
		$gene[$i]= $gene[$i];
		$start[$i]= $start[$i];
		$end[$i]= $end[$i];
		print de_dup "$OP[$i]\t$strand[$i]\t$gene[$i]\t$start[$i]\t$end[$i]\n";
	}
}
no warnings 'uninitialized';	
open(fh2,"de_dup_operon_struc_door");
#print "file opened for reading \n";
my @read2=<fh2>;
my $count_read2=0;
foreach (@read2) {
	chomp;
	$count_read2++;
}
@OP=0; @start=0; @end=0; @gene=0; @strand=0;
my @values1=0;
for(my $i=0; $i<$count_read2; $i++){
	@values1= split("\t", $read2[$i]);
	$OP[$i]= $values1[0];
	$start[$i]= $values1[3];
	$end[$i]= $values1[4];
	$gene[$i]= $values1[2];
	$strand[$i]= $values1[1];
}
open(de_dup_next, '>operon_map_door');

$b=0;
for (my $i=0; $i < $count_read2; $i++){
	if($OP[$i+1]==$OP[$i]){
		$b=1;
	}
	else{
		print de_dup_next "$OP[$i]\t$gene[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\n";
	}
}
	
system_commands1();

sub system_commands1{

	system('rm de_dup_operon_struc_door operon_struc_door');
	}
	

#!/usr/bin/perl
#############################################
#SLOFE calculates the ratio for an operon using the deltaG (free-energy) of the harbored stem-loops. 
#The calculated ratio is simply the ratio of the deltaG of stem-loops. 
#This ratio represents the stoichiometry of the encoded complex by that operon.
# 
#Author: Yogendra Bhaskar
#
#License: GPL 3.0
#
#Usage: perl SLOFE.pl -h
#############################################

use strict;
use warnings;

if(!@ARGV){
	print "No FilePath Specified!\n";
	exit;
}

if($ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
	help();
	exit;
}

(@ARGV == 4) || die "SLOFE runs with four argument.\nHelp: perl SLOFE -h\nUsage: perl SLOFE <genome_fasta_file> <descr> <annotation_file> <operon_map>\n";

sub help { print "
                                SLOFE v1.0\n
Usage: ./SLOFE <genome_fasta_file> <descr> <annotation_file> <operon_map>\n
SLOFE takes genome file in fasta format, descr file (readme), gene annotation and operon annotation.
1. Genome file must have the extension '.fasta'.
2. Descriptor (descr) file is used by RNAmotif to predict the desired motifs, which is used for RNA secondary stucture prediction. 
   Descriptor file must have the extension '.descr'.
3. Gene annotation can be provided in the same format as given in sample_input folder or it can be downloaded from the genbank. 
   The downloaded file must be prepared in the accepeted format using the script 'annotation_prep' (readme).
4. Operon annotation/map can be downloaded either from Door database or ProOpDb database. 
   Then it should be prepared in the accepted fromat using scripts 'operon-Door' or 'operon-ProOpDb' (readme).

Example: perl SLOFE sample_input/nc.fasta sample_input/H10-mod.descr sample_input/gene_annotation sample_input/operon_map
\n";}


my $fasta = $ARGV[0];
my $descr = $ARGV[1];
my $GB = $ARGV[2];
my $Operon= $ARGV[3];


die "Input 1 must have extension .fasta" if ($fasta !~ /\.fasta/);
die "Input 2 must have extension .descr" if ($descr !~ /\.descr/);

my $tmp="tmp_slofe";
if (-d $tmp){
	print "Temporary files will be written in \"tmp_slofe/\" directory.\n";# directory called tmp_slofe exists
}
else{
	mkdir $tmp;
	print "Temporary files will be written in \"tmp_slofe/\" directory.\n";
}

system_commands1($fasta, $descr);

sub system_commands1{
	print "Motif Prediction\n";
	system('./rnamotif -descr ' .$descr.' '. $fasta.' > tmp_slofe/rnamotif.out');
	
	if(-z "tmp_slofe/rnamotif.out"){ die "rnamotif.out doesn't exist or has zero size";}	

#	print "Redundancy removal using rmprune\n";
	system('./rmprune tmp_slofe/rnamotif.out > tmp_slofe/rmprune1.out');
	system('./rmprune tmp_slofe/rmprune1.out > tmp_slofe/rmprune2.out');
	system('./rmprune tmp_slofe/rmprune2.out > tmp_slofe/rmprune.out');

	if(-z "tmp_slofe/rmprune.out"){ die "rmprune.out doesn't exist or has zero size";}

#	print "Sorting\n";
	system("sed -e '/^>/d' -e '/#/d' tmp_slofe/rmprune.out | sort -k4 -g -u > tmp_slofe/rmprune.out_uniq");

	if(-z "tmp_slofe/rmprune.out_uniq"){ die "rmprune.out_uniq doesn't exist or has zero size";}

#	print "Extracting sequences using awk\n";
	system('awk \'{FS=" "}{ $1=$2=$3=$4=$5=""; print $0}\' tmp_slofe/rmprune.out_uniq|tr -d \' \' > tmp_slofe/sequence.out');

	my $rnamotif_uniq= `wc -l tmp_slofe/sequence.out`;
	print "Number of predicted motifs: $rnamotif_uniq";

	if(-z "tmp_slofe/sequence.out"){ die "sequence.out doesn't exist or has zero size";}

	print "RNA structure calculation using RNAfold\n";
	system('./RNAfold -d2 --noLP < tmp_slofe/sequence.out > tmp_slofe/rnafold.out');

	if(-z "tmp_slofe/rnafold.out"){ die "rnafold.out doesn't exist or has zero size";}

#	print "Modifying RNAfold output\n";
	system('sed -e \'/^[.,(]/s/$/@/\' tmp_slofe/rnafold.out |tr "\n" "\t"|tr "@" "\n" > tmp_slofe/rnafold.out_modi');

	if(-z "tmp_slofe/rnafold.out_modi"){ die "rnafold.out_modi doesn't exist or has zero size";}

#	print "Combining RNAmotif and RNAfold output\n";
	system('sed \'s/^\\t//g\' tmp_slofe/rnafold.out_modi |sed \'1i\\\\\' |awk \'{FS="\\t"}{  print $1}\'|sed \'1d\' > tmp_slofe/seq_from_rnafold');
	system('awk \'{FS=" "}{  print $2}\' tmp_slofe/rnafold.out_modi > tmp_slofe/fold_from_rnafold');
	system('awk \'{FS="("}{print $NF}\' tmp_slofe/rnafold.out_modi|sed \'s/)//g\' > tmp_slofe/eng_from_rnafold');
	system('paste tmp_slofe/rmprune.out_uniq tmp_slofe/seq_from_rnafold tmp_slofe/fold_from_rnafold tmp_slofe/eng_from_rnafold |awk \'{FS=" "}{print $1"\t"$3"\t"$4"\t"$5"\t"$6$7$8"\t"$9"\t"$10"\t"$11}\'> tmp_slofe/rnaM_rnaF_merged');
	print "RNAfold prediction completed\n";
	if(-z "tmp_slofe/rnaM_rnaF_merged"){ die "rnaM_rnaF_merged doesn't exist or has zero size";}

#	print "deltaG cutoff of -5\n";
	system('awk \'{FS="\t"}{if($8<=-5) print $0}\' tmp_slofe/rnaM_rnaF_merged > tmp_slofe/cutoff5.out');
	system('sed \'1i\\\\\' tmp_slofe/cutoff5.out |awk \'{FS=" "}{if($2==0){ print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($3+$4-1)} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"($3-$4+1)}}\'|sed \'1d\'|sed \'$d\' > tmp_slofe/cutoff5_new.csv');

	if(-z "tmp_slofe/cutoff5_new.csv"){ die "cutoff5_new.csv doesn't exist or has zero size";}	

#	system('grep -v \').*(\' tmp_slofe/cutoff5_new.csv |awk \'{FS="\t"}{if($2==0) {print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print $1"\t"$2"\t"$9"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}\'|sort -k3 -g > tmp_slofe/cutoff5_no2braces_sorted');
	system('grep -v \').*(\' tmp_slofe/cutoff5_new.csv |awk \'{FS="\t"}{if($2==0) {print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print $1"\t"$2"\t"$9"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}\'|sort -k3 -g > tmp_slofe/cutoff5_stem-loops');

	my $cutoff5= `wc -l tmp_slofe/cutoff5_stem-loops`;
	print "Number of motifs after -5 deltaG cutoff: $cutoff5";
}

if(-z "tmp_slofe/cutoff5_stem-loops"){ die "cutoff5_stem-loops doesn't exist or has zero size";}


#print "Removing 3 plus overlapped sequences\n";

open(fh1,"tmp_slofe/cutoff5_stem-loops");
#print "file opened for reading \n";
my @read1=<fh1>;
my $count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

my @name2; my @strand2; my @start2; my @end2; my @length2; my @rm_seq2; my @rf_seq2; my @fold2; my @energy2;

for(my $i=0; $i<$count_read1; $i++){
	my @values1= split("\t", $read1[$i]);
	$name2[$i]= $values1[0];
	$strand2[$i]= $values1[1];
	$start2[$i]= $values1[2];
	$end2[$i]= $values1[3];
	$length2[$i]= $values1[4];
	$rm_seq2[$i]= $values1[5];
	$rf_seq2[$i]= $values1[6];
	$fold2[$i]= $values1[7];
	$energy2[$i]= $values1[8];
}

open(Removed, '>tmp_slofe/3plus_removed');
my $plus3_removed=0;
for(my $i=0; $i< $count_read1; $i++) {
	if($start2[$i] >= $start2[$i-1] && $end2[$i] <= ($end2[$i-1] +2)){
		$plus3_removed=1;
#		print "$name[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\tF\n";
	}
	else{

		print Removed "$name2[$i]\t$strand2[$i]\t$start2[$i]\t$end2[$i]\t$length2[$i]\t$rm_seq2[$i]\t$rf_seq2[$i]\t$fold2[$i]\t$energy2[$i]\n";	
	}
}
close(fh1);
close(Removed);

open(fh2,"tmp_slofe/3plus_removed");


#print "file opened for reading \n";
my @read2=<fh2>;
my $count_read2=0;
foreach (@read2) {
	chomp;
	$count_read2++;
}

my (@name1, @strand1, @start1, @end1, @length1, @rm_seq1, @rf_seq1, @fold1, @energy1);

for(my $i=0; $i<$count_read2; $i++){
	my @values2= split("\t", $read2[$i]);
	$name1[$i]= $values2[0];
	$strand1[$i]= $values2[1];
	$start1[$i]= $values2[2];
	$end1[$i]= $values2[3];
	$length1[$i]= $values2[4];
	$rm_seq1[$i]= $values2[5];
	$rf_seq1[$i]= $values2[6];
	$fold1[$i]= $values2[7];
	$energy1[$i]= $values2[8];
}

open(Removed1, '>tmp_slofe/non-redundant_stem-loops');
my $plus3_removed1=0;
for(my $i=0; $i< $count_read2; $i++) {
	if($start1[$i] >= $start1[$i-1] && $end1[$i] <= ($end1[$i-1] +2)){
		$plus3_removed1=1;
#		print "$name1[$i]\t$strand1[$i]\t$start1[$i]\t$end1[$i]\t$length1[$i]\t$rm_seq1[$i]\t$rf_seq1[$i]\t$fold1[$i]\t$energy1[$i]\tF\n";
	}
	else{
		print Removed1 "$name1[$i]\t$strand1[$i]\t$start1[$i]\t$end1[$i]\t$length1[$i]\t$rm_seq1[$i]\t$rf_seq1[$i]\t$fold1[$i]\t$energy1[$i]\n";
	}
}
close(fh2);
close(Removed1);
if(-z "tmp_slofe/non-redundant_stem-loops"){ die "non-redundant_stem-loops doesn't exist or has zero size";}

my $uniq_stemloops= `wc -l tmp_slofe/non-redundant_stem-loops`;
	print "Number of predicted stem-loops: $uniq_stemloops";

system ('awk \'{FS="\t"}{if($2==0) print $0}\' tmp_slofe/non-redundant_stem-loops > tmp_slofe/for_G_map_0_improved');
system ('awk \'{FS="\t"}{if($2==1) print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}\' tmp_slofe/non-redundant_stem-loops |sort -k4 -g > tmp_slofe/sorted_for_G_map_1_improved');

print "Genome mapping of postive strand stem-loops\n";

open(fh3,"tmp_slofe/for_G_map_0_improved");
#print "file opened for reading \n";
my @read3=<fh3>;
my $count_read3=0;
foreach (@read3) {
	chomp;
	$count_read3++;
}

my @name; my @strand; my @start; my @end; my @length; my @rm_seq; my @rf_seq; my @fold; my @energy;

for(my $i=0; $i<$count_read3; $i++){
	my @values3= split("\t", $read3[$i]);
	$name[$i]= $values3[0];
	$strand[$i]= $values3[1];
	$start[$i]= $values3[2];	
	$end[$i]= $values3[3];
	$length[$i]= $values3[4];
	$rm_seq[$i]= $values3[5];
	$rf_seq[$i]= $values3[6];
	$fold[$i]= $values3[7];
	$energy[$i]= $values3[8];
} 

open(fh4, $GB);
#print "file opened for reading  \n";
my @read4=<fh4>;
my $count_read4 = 0;
foreach (@read4) {
	chomp;
	$count_read4++;
}
#print "number of ids: $count_read\n";

my (@gene_start, @gene_end, @gene_strand, @gene_new);

for (my $i=0; $i<$count_read4; $i++) { 
	my @values4 = split("\t",$read4[$i]);	
	$gene_start[$i]= $values4[0];
	$gene_end[$i]= $values4[1];
	$gene_strand[$i]= $values4[3];
	$gene_new[$i]= $values4[2];
}

open(Gmap0, '>tmp_slofe/g_mapped_0');
no warnings 'uninitialized';
for(my $i=0; $i< $count_read3; $i++) {
	for(my $j=0; $j< $count_read4; $j++) {
		if($start[$i] >= $gene_start[$j] && $end[$i] <= $gene_end[$j]) {
			print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\tinside->$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
		elsif ($start[$i] >= $gene_start[$j] && $start[$i] <= $gene_end[$j] &&$end[$i] >= $gene_end[$j] &&$end[$i] <= $gene_start[$j+1]) {
			if($gene_strand[$j] eq "-"){
				print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_5'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j]/$gene_strand[$j]/_and_$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
			else{
				print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_3'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j]/$gene_strand[$j]/_and_$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
		}
		elsif ($start[$i] > $gene_end[$j] && $end[$i] < $gene_start[$j+1]) {
			print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\tintergenic->$gene_new[$j]/$gene_strand[$j]/<->$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
		elsif ($start[$i] <= $gene_start[$j] && $end[$i] >= $gene_start[$j] && $end[$i] <= $gene_end[$j]) {
			if($gene_strand[$j] eq "-"){
				print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_3'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j-1]/$gene_strand[$j-1]/_and_$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
			else{
				print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_5'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j-1]/$gene_strand[$j-1]/_and_$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
		}
		elsif($start[$i] >= $gene_start[$j] && $start[$i] <= $gene_end[$j] && $end[$i] >= $gene_start[$j+1] && $end[$i] <= $gene_end[$j+1]) {
			print Gmap0 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_with->$gene_new[$j]*$gene_strand[$j]*_and_$gene_new[$j+1]*$gene_strand[$j+1]*\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
	}
}

close(fh3);
#close(fh4);
close(Gmap0);

print "Genome mapping of negative strand stem-loops\n";

open(fh5,"tmp_slofe/sorted_for_G_map_1_improved");
#print "file opened for reading \n";
my @read5=<fh5>;
my $count_read5=0;
foreach (@read5) {
	chomp;
	$count_read5++;
}

@name=0; @strand=0; @start=0;  @end=0;  @length=0;  @rm_seq=0;  @rf_seq=0;  @fold=0;  @energy=0;

for(my $i=0; $i<$count_read5; $i++){
	my @values5= split("\t", $read5[$i]);
	$name[$i]= $values5[0];
	$strand[$i]= $values5[1];
	$start[$i]= $values5[2];	
	$end[$i]= $values5[3];
	$length[$i]= $values5[4];
	$rm_seq[$i]= $values5[5];
	$rf_seq[$i]= $values5[6];
	$fold[$i]= $values5[7];
	$energy[$i]= $values5[8];
}

#open(fh6, "annotation_frm_gb");
#print "file opened for reading  \n";
#@read2=<fh2>;
#$count_read2 = 0;
#foreach (@read2) {
#	chomp;
#	$count_read2++;
#}
#print "number of ids: $count_read\n";

#for ($i=0; $i<$count_read2; $i++) { 
#	@values = split("\t",$read2[$i]);	
#	$gene_start[$i]= $values[0];
#	$gene_end[$i]= $values[1];
#	$gene_strand[$i]= $values[3];
#	$gene_new[$i]= $values[2];
#}

open(Gmap1, '>tmp_slofe/g_mapped_1');

for(my $i=0; $i< $count_read5; $i++) {
	for(my $j=0; $j< $count_read4; $j++) {
		if($start[$i] <= $gene_end[$j] && $end[$i] >= $gene_start[$j]) {
			print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\tinside->$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
		elsif ($start[$i] >= $gene_end[$j] && $end[$i] <= $gene_end[$j] &&$end[$i] >= $gene_start[$j] &&$start[$i] <= $gene_start[$j+1]) {
			if($gene_strand[$j] eq "-"){
				print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_5'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j]/$gene_strand[$j]/_and_$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
			else{
				print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_3'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j]/$gene_strand[$j]/_and_$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
		}
		elsif ($end[$i] > $gene_end[$j] && $start[$i] < $gene_start[$j+1]) {
			print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\tintergenic->$gene_new[$j]/$gene_strand[$j]/<->$gene_new[$j+1]/$gene_strand[$j+1]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
		elsif ($end[$i] <= $gene_start[$j] && $start[$i] >= $gene_start[$j] && $start[$i] <= $gene_end[$j]) {
			if($gene_strand[$j] eq "-"){
				print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_3'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j-1]/$gene_strand[$j-1]/_and_$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
			else{
				print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_on_5'_with->$gene_new[$j]/$gene_strand[$j]/_btw_$gene_new[$j-1]/$gene_strand[$j-1]/_and_$gene_new[$j]/$gene_strand[$j]/\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
			}
		}
		elsif($end[$i] >= $gene_start[$j] && $end[$i] <= $gene_end[$j] && $start[$i] >= $gene_start[$j+1] && $start[$i] <= $gene_end[$j+1]) {
			print Gmap1 "$gene_start[$j]\t$gene_end[$j]\t$gene_strand[$j]\t$gene_new[$j]\t$strand[$i]\t$start[$i]\t$end[$i]\toverlaped_with->$gene_new[$j]*$gene_strand[$j]*_and_$gene_new[$j+1]*$gene_strand[$j+1]*\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
		}
	}
}

close(fh4);
close(fh5);
close(Gmap1);

system_commands2();

sub system_commands2{
#	print "length cut according to fold\n";
	system('cat tmp_slofe/g_mapped_0  tmp_slofe/g_mapped_1 |sort -k6 -g > tmp_slofe/g_mapped_both');
	system('cut -f12 tmp_slofe/g_mapped_both | cut -d\\( -f1|awk \'{print length}\' > tmp_slofe/left_cut_length');	
	system('cut -f12 tmp_slofe/g_mapped_both | rev| cut -d\\) -f1| awk \'{print length}\' > tmp_slofe/right_cut_length ');
	system('cut -f5,6 tmp_slofe/g_mapped_both > tmp_slofe/tmp_start ');
	system('cut -f5,7 tmp_slofe/g_mapped_both > tmp_slofe/tmp_end');
	system('paste tmp_slofe/tmp_start tmp_slofe/left_cut_length | awk \'{FS="\t"}{if($1==0) {print $2+$3} else if($1==1){print $2-$3}}\' > tmp_slofe/start');
	system('paste tmp_slofe/tmp_end tmp_slofe/right_cut_length | awk \'{FS="\t"}{if($1==0) {print $2-$3} else if($1==1){print $2+$3}}\' > tmp_slofe/end');
	system('cut -f10 tmp_slofe/g_mapped_both|sed \'s/ //g\' > tmp_slofe/seq_from_g_mapped_both');
	system('paste tmp_slofe/seq_from_g_mapped_both tmp_slofe/left_cut_length |awk \'{FS="\t"}{print substr($1,($2+1),length($1)-$2)}\' > tmp_slofe/left_already_cut_seq');
	system('paste tmp_slofe/left_already_cut_seq tmp_slofe/right_cut_length |awk \'{FS="\t"}{print substr($1,1,length($1)-$2)}\' > tmp_slofe/new_stem_seq');

#	print "new start, end and stem_seq\n";
	system('paste tmp_slofe/g_mapped_both tmp_slofe/start tmp_slofe/end tmp_slofe/new_stem_seq|sed \'1i\\\\\'|awk \'{FS="\t"}{if($5==0) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$14"\t"$15"\t"$8"\t"length($16)"\t"$16"\t"$11"\t"$12"\t"$13} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$15"\t"$14"\t"$8"\t"length($16)"\t"$16"\t"$11"\t"$12"\t"$13}}\'|sed \'1d\'|sort -k6 -g -u > tmp_slofe/g_mapped_modi');

#	print "Right stem length\n";
	system('cut -f12 tmp_slofe/g_mapped_modi|sed -e \'s/(//g\' -e \'s/\\.//g\'|awk \'{print length}\' > tmp_slofe/right_stem_length');

	print "Calculating eng*stem/stem-loop\n";
	system('paste tmp_slofe/g_mapped_modi tmp_slofe/right_stem_length |awk \'{FS="\t"}{print $0"\t"$13/$9"\t"$13*$14"\t"($13*$14)/$9}\' > tmp_slofe/g_mapped_modi1 ');

#	print "searched the intergenic||overlaped and add 1 in the end otherwise 0. Sorted based on the (energy*stem_length)/stem-loop_length\n";
	system('awk \'{FS="\t"}{if($8~/intergenic|overlaped/) {print $0"\t1"} else {print $0"\t0"}}\' tmp_slofe/g_mapped_modi1 |sort -k17 -g  > tmp_slofe/for_screening_calc');
	
#	print "did sum of the column and calculated the percent. Then selected the over 60% values\n";
	system('awk \'{FS="\t"}{print $0"\t"(sum+=$18)"\t"NR}\' tmp_slofe/for_screening_calc|awk \'{FS="\t"}{print $0"\t"($19/$20)}\'|awk \'{FS="\t"}{if($21>0.59) print $0}\' > tmp_slofe/final_stem-loops');

	system('sort -k6 -g -u tmp_slofe/final_stem-loops > tmp_slofe/final_stem-loops_sorted');
	my $final_stemloops= `wc -l tmp_slofe/final_stem-loops_sorted`;
	print "Number of final stem-loops after threshold based on staility factor: $final_stemloops";

}

#print "filtered final stem-loops\n";

open(fh6,"tmp_slofe/final_stem-loops_sorted");
#print "file opened for reading \n";
my @read6=<fh6>;
my $count_read6=0;
foreach (@read6) {
	chomp;
	$count_read6++;
}

my @g_start=0; my @g_end=0; my @g_strand=0; my @gene=0; @strand=0;  @start=0;  @end=0; my @region=0;  @length=0;  @rm_seq=0;  @rf_seq=0;  @fold=0;  @energy=0;

for(my $i=0; $i<$count_read6; $i++){
	my @values6= split("\t", $read6[$i]);
	$g_start[$i]= $values6[0];
	$g_end[$i]= $values6[1];
	$g_strand[$i]= $values6[2];
	$gene[$i]=$values6[3];
	$strand[$i]= $values6[4];
	$start[$i]= $values6[5];
	$end[$i]= $values6[6];
	$region[$i]= $values6[7];
	$length[$i]= $values6[8];
	$rm_seq[$i]= $values6[9];
	$rf_seq[$i]= $values6[10];
	$fold[$i]= $values6[11];
	$energy[$i]= $values6[12];
}

open(filter1, '>tmp_slofe/final_stem-loops_filtered');

print filter1 "$g_start[0]\t$g_end[0]\t$g_strand[0]\t$gene[0]\t$strand[0]\t$start[0]\t$end[0]\t$region[0]\t$length[0]\t$rm_seq[0]\t$rf_seq[0]\t$fold[0]\t$energy[0]\n";

for(my $i=0; $i< $count_read6; $i++) {
	if($start[$i] >= $start[$i-1] && $end[$i] >= $end[$i-1]){
		print filter1 "$g_start[$i]\t$g_end[$i]\t$g_strand[$i]\t$gene[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$region[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$fold[$i]\t$energy[$i]\n";
	}
}

close(fh6);
close(filter1);


print "Operon mapping\n";

open(fh7,"tmp_slofe/final_stem-loops_filtered");
#print "file opened for reading \n";
my @read7=<fh7>;
my $count_read7=0;
foreach (@read7) {
	chomp;
	$count_read7++;
}

@g_start=0;  @g_end=0;  @g_strand=0;  @gene=0;  @strand=0;  @start=0;  @end=0;  @region=0;  @length=0;  @rm_seq=0;  @rf_seq=0;  @fold=0;  @energy=0; my @folding=0;
my @values1=0;

for(my $i=0; $i<$count_read7; $i++){
	@values1= split("\t", $read7[$i]);	
	$g_start[$i]= $values1[0];
	$g_end[$i]= $values1[1];
	$g_strand[$i]= $values1[2];
	$gene_new[$i]= $values1[3];	
	$strand[$i]= $values1[4];
	$start[$i]= $values1[5];	
	$end[$i]= $values1[6];
	$region[$i]= $values1[7];
	$length[$i]= $values1[8];
	$rm_seq[$i]= $values1[9];
	$rf_seq[$i]= $values1[10];
	$energy[$i]= $values1[12];
	$folding[$i]= $values1[11];
}


open(fh8, $Operon);
#print "file opened for reading  \n";
my @read8=<fh8>;
my $count_read8 = 0;
foreach (@read8) {
	chomp;
	$count_read8++;
}
#print "number of ids: $count_read\n";

my @operon=0; my @gene_new1=0; my @gene_old1=0;  @gene_start=0;  @gene_end=0;  @gene_strand=0;
my @values=0;

for (my $i=0; $i<$count_read8; $i++) { 
	@values = split("\t",$read8[$i]);
	$operon[$i]= $values[0];	
	$gene_new1[$i]= $values[1];
	$gene_old1[$i]= $values[2];	
	$gene_start[$i]= $values[4];
	$gene_end[$i]= $values[5];
	$gene_strand[$i]= $values[3];
}

open(OP, '>tmp_slofe/final_stem-loops_operon_mapped');

for(my $i=0; $i< $count_read7; $i++) {
	for(my $j=0; $j< $count_read8; $j++) {
			if($start[$i] >= $gene_start[$j] && $end[$i] <= $gene_end[$j]){
				print OP "$g_start[$i]\t$g_end[$i]\t$g_strand[$i]\t$gene_new[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$region[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$folding[$i]\t$energy[$i]\t$operon[$j]\tinsidOP->$operon[$j]\n";
			}
			elsif($start[$i] > $gene_end[$j] && $end[$i] < $gene_start[$j+1]) {
				print OP "$g_start[$i]\t$g_end[$i]\t$g_strand[$i]\t$gene_new[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$region[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$folding[$i]\t$energy[$i]\t$operon[$j]\tinteroperonic->$operon[$j]<->$operon[$j+1]\n";
			}
			elsif ($start[$i] <= $gene_start[$j] && $end[$i] >= $gene_start[$j] && $end[$i] <= $gene_end[$j]) {
				print OP "$g_start[$i]\t$g_end[$i]\t$g_strand[$i]\t$gene_new[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$region[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$folding[$i]\t$energy[$i]\t$operon[$j]\toverlapOP->$operon[$j]_btw_$operon[$j-1]<->$operon[$j]\n";
			}
			elsif ($start[$i] >= $gene_start[$j] && $start[$i] <= $gene_end[$j] && $end[$i] >= $gene_end[$j] && $end[$i] <= $gene_start[$j+1]) {
				print OP "$g_start[$i]\t$g_end[$i]\t$g_strand[$i]\t$gene_new[$i]\t$strand[$i]\t$start[$i]\t$end[$i]\t$region[$i]\t$length[$i]\t$rm_seq[$i]\t$rf_seq[$i]\t$folding[$i]\t$energy[$i]\t$operon[$j]\toverlapOP->$operon[$j]_btw_$operon[$j]<->$operon[$j+1]\n";
			}
	
	}
}

close(fh7);
close(fh8);
close(OP);

system_commands3();

sub system_commands3{
#	print "T extended sequence calculation\n";
	system('cut -f5,6,7 tmp_slofe/final_stem-loops_operon_mapped > tmp_slofe/coords');
	system('perl extract_seq.pl '.$fasta.'  tmp_slofe/coords > tmp_slofe/T_ext');
	system('cut -f12 tmp_slofe/final_stem-loops_operon_mapped |sed -e \'s/\\.*$//g\'|sed -e \'s/^\\.*//g\' > tmp_slofe/fold_no-outside_dots');
	system('sed \'1i\\\\\' tmp_slofe/fold_no-outside_dots |awk \'{FS="("} {print $NF}\'|rev|awk \'{FS=")"}{print length"\t"length($NF)"\t"length-length($NF)}\'|sed \'1d\'|cut -f3 > tmp_slofe/right_stem_actual_length');
	system('cut -f10 tmp_slofe/final_stem-loops_operon_mapped > tmp_slofe/stem-loop_seq');
	system('paste tmp_slofe/stem-loop_seq tmp_slofe/right_stem_actual_length |awk \'{print substr($1, length($1)-($2-1),$2)}\' > tmp_slofe/right_stem_seq');
	system('cut -f10 tmp_slofe/final_stem-loops_operon_mapped |tr t u > tmp_slofe/new_rf_seq');
	system('cut -f5 tmp_slofe/T_ext|tr ATGC atgc > tmp_slofe/extended10_seq');
	system('paste tmp_slofe/final_stem-loops_operon_mapped tmp_slofe/fold_no-outside_dots tmp_slofe/right_stem_actual_length tmp_slofe/extended10_seq tmp_slofe/right_stem_seq tmp_slofe/new_rf_seq |cut -f-10,13-|awk \'{FS="\t"}{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$18"\t"$11"\t"$14"\t"$15"\t"$16"\t"$17$16"\t"$12"\t"$13}\' > tmp_slofe/final_stem-loops_improved');
	print "T content and Poly U\n";	
	system('cut -f16 tmp_slofe/final_stem-loops_improved > tmp_slofe/extended_right_stem_seq');
	system('perl maxTcontent.pl tmp_slofe/extended_right_stem_seq |awk \'{FS="\t"}{if($3=="" && $4=="" && $5==""){print $0"\tNo-PolyU"} else if($4=="" && $5=="") {print $0"\tPoly3U"} else if($3!="" && $4!="" && $5==""){print $0"\tPoly4U"} else if($5!=""){print $0"\tPoly5U"}}\' > tmp_slofe/T_cont_polyT');
	system('cut -f2,6 tmp_slofe/T_cont_polyT > tmp_slofe/Poly_T');
#	system('paste tmp_slofe/final_stem-loops_improved tmp_slofe/Poly_T | grep -v \').*(\' > final_stem-loops_improved_T');
	system('paste tmp_slofe/final_stem-loops_improved tmp_slofe/Poly_T | grep -v \').*(\' > stable_stem-loops');
	my $stable_stemloops= `wc -l stable_stem-loops`;
	print "Number of stable stem-loops: $stable_stemloops";
	system('grep insidOP stable_stem-loops |grep -v inside > intergenic_intraOperonic_stem-loops');
	print "Stable stem-loops prediction finished\n";
	print "Ratio Calculation.......\n";	
	system('perl ratio_calc.pl  '.$Operon.' '.$GB.'');
}

__END__

#!/usr/bin/perl
#############################################
#This program is a part of the SLOFE and takes inputs from it.
#This program indentifies the SRPS operons and calculates the ratio using the predicted stable stem-loops.
#
#Author: Yogendra Bhaskar
#
#License: GPL 3.0
#
#############################################

use strict;
use warnings;
use List::Util qw( min max ); #which is released with the Perl distribution v5.7.3

my $Operon= $ARGV[0];
my $GB = $ARGV[1];

system_commands($Operon, $GB);

open(fh2,"tmp_slofe/srps_operon_genes_new");
#print "file opened for reading \n";
my @read2=<fh2>;
my $count_read2=0;
foreach (@read2) {
	chomp;
	$count_read2++;
}

my @OP2; my @strand2; my @op_count2; my @gene2; my @g_start2; my @g_end2;

for(my $i=0; $i<$count_read2; $i++){
	my @values2= split("\t", $read2[$i]);	
	$OP2[$i] = $values2[0];
	  $strand2[$i]=  $values2[1];
	  $op_count2[$i]=  $values2[2];
	  $gene2[$i]=  $values2[3];
	  $g_start2[$i]=  $values2[4];
	  $g_end2[$i]=  $values2[5];
}

close(fh2);

open(fh3,"tmp_slofe/final_stem-loops_improved_T_intergenic");
#print "file opened for reading \n";
my @read3=<fh3>;
my $count_read3=0;
foreach (@read3) {
	chomp;
	$count_read3++;
}

my @g_start3; my @g_end3; my @g_strand3; my @gene_new3; my @strand3; my @start3; my @end3; my @region3; my @length3; my @rm_seq3; my @rf_seq3; my @energy3; my @folding3; my @OP3; my @OP_region3; my @U_c3; my @PolyU3;

for(my $i=0; $i<$count_read3; $i++){
	my @values3= split("\t", $read3[$i]);	
	$g_start3[$i]= $values3[0];
	$g_end3[$i]= $values3[1];
	$g_strand3[$i]= $values3[2];
	$gene_new3[$i]= $values3[3];	
	$strand3[$i]= $values3[4];
	$start3[$i]= $values3[5];	
	$end3[$i]= $values3[6];
	$region3[$i]= $values3[7];
	$length3[$i]= $values3[8];
	$rm_seq3[$i]= $values3[9];
	$rf_seq3[$i]= $values3[10];
	$energy3[$i]= $values3[11];
	$folding3[$i]= $values3[12];
	$OP3[$i]= $values3[16];
	$OP_region3[$i]= $values3[17];
	$U_c3[$i]= $values3[18];
	$PolyU3[$i]= $values3[19];	
}
close(fh3);

no warnings 'uninitialized';
my @eng_new;
my @reg_new;
for(my $i=0; $i<$count_read2; $i++){
	for(my $j=0; $j<$count_read3; $j++){
		if($strand2[$i] eq '-'){
			if($OP2[$i] == $OP3[$j] && $gene2[$i] eq $gene_new3[$j] && $OP2[$i] == $OP2[$i+1]){
				$eng_new[$i+1]=$energy3[$j];
				$reg_new[$i+1]=$region3[$j];
#				print "$OP2[$i]\t$strand2[$i]\t$op_count2[$i]\t$gene2[$i]\t$g_start2[$i]\t$g_end2[$i]\t$eng_new[$i]\tmid\n";
			}
			elsif($OP2[$i] != $OP2[$i-1] && $OP3[$j] == ($OP2[$i] -1) && $OP_region3[$j] =~ m/$OP2[$i]/){
				$eng_new[$i]=$energy3[$j];
				$reg_new[$i]=$region3[$j];
#				print "$OP2[$i]\t$strand2[$i]\t$op_count2[$i]\t$gene2[$i]\t$g_start2[$i]\t$g_end2[$i]\t$eng_new[$i]\tend\n";
			}
		}
		elsif($strand2[$i] eq '+'){
			if($OP2[$i] == $OP3[$j] && $gene2[$i] eq $gene_new3[$j]){
				$eng_new[$i]=$energy3[$j];
				$reg_new[$i]=$region3[$j];				
			}
		}
	}
}

open (wrfile1, '>tmp_slofe/with_eng_redef');
for(my $i=0; $i<$count_read2; $i++){
	if($reg_new[$i] eq "" && $eng_new[$i] eq ""){
		print wrfile1 "$OP2[$i]\t$strand2[$i]\t$op_count2[$i]\t$gene2[$i]\t$g_start2[$i]\t$g_end2[$i]\t0\t0\n";
	}
	else{
		print wrfile1 "$OP2[$i]\t$strand2[$i]\t$op_count2[$i]\t$gene2[$i]\t$g_start2[$i]\t$g_end2[$i]\t$reg_new[$i]\t$eng_new[$i]\n";
	}
}
close(wrfile1);

open(fh1,"tmp_slofe/with_eng_redef");
#print "file opened for reading \n";
my @read1=<fh1>;
my $count_read1=0;
foreach (@read1) {
	chomp;
	$count_read1++;
}

my @values1;
my (@OP, @strand,@op_count,@gene,@g_start,@g_end,@region,@eng);

for(my $i=0; $i<$count_read1; $i++){
	@values1= split("\t", $read1[$i]);	
	$OP[$i] = $values1[0];
	  $strand[$i]=  $values1[1];
	  $op_count[$i]=  $values1[2];
	  $gene[$i]=  $values1[3];
	  $g_start[$i]=  $values1[4];
	  $g_end[$i]=  $values1[5];
	  $region[$i]=  $values1[6];
	  $eng[$i]=  $values1[7];
#	  print "$OP[$i]\t$strand[$i]\t$op_count[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\n";	
}

#Ratio calculation

my $a=1 ;
until($a>int((0.8*max(@op_count))+0.99)){
	for(my $i=0; $i<$count_read1; $i++){
		if($strand[$i] eq '+' && $OP[$i]==$OP[$i+1]){
			if($eng[$i]==0){
#				print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i+1]\tnewP\n";
				$eng[$i]=$eng[$i+1];
			}
			else{
#				print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\torgP\n";	
				$eng[$i]=$eng[$i];
			}
		}
		elsif($strand[$i] eq '-' && $OP[$i]==$OP[$i-1]){
			if($eng[$i]==0){
#				print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i-1]\tnewN\n";
				$eng[$i]=$eng[$i-1];
			}
			else{
#				print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\torgN\n";
				$eng[$i]=$eng[$i];
			}
		}
		else{
#			print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\tNEOP\n";
			$eng[$i]=$eng[$i];
		}
	}
	$a=$a+1;	
}

my @eng2=0;

for(my $i=0; $i<$count_read1; $i++){
	if($strand[$i] eq '+' && $OP[$i]==$OP[$i+1]){
		if($strand[$i] eq '+' && $OP[$i]!=$OP[$i-1]){
			$eng2[$i]=$eng[$i];
		}
#		$eng2[$i]=$eng[$i];
		$eng2[$i+1]=$eng2[$i];
#		$eng2[$i]= $eng2[$i+1];
	}
#	elsif($strand[$i] eq '-' && $OP[$i]==$OP[$i-1]){
#		$eng2[$i]=$eng[$i];
#		$eng2[$i-1]=$eng[$i];
#	}
	if($OP[$i]!=$OP[$i-1] && $OP[$i]!=$OP[$i+1]){
		$eng[$i]=1;
		$eng2[$i]=1;
	}
}

for(my $i=$count_read1; $i>0; $i--){
	if($strand[$i] eq '-' && $OP[$i]==$OP[$i-1]){
		if($strand[$i] eq '-' && $OP[$i]!=$OP[$i+1]){
			$eng2[$i]=$eng[$i];
		}
		$eng2[$i-1]=$eng2[$i];
	}	
}

# Ratio normalization 
my @ratio=0;
for(my $i=0; $i<$count_read1; $i++){
	$ratio[$i]= sprintf("%.2f", abs($eng[$i]/$eng2[$i]));
#	print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\t$eng2[$i]\t$ratio[$i]\n";
}

#Ratio representation in this form: 1:0:0:0

for(my $i=0; $i<$count_read1; $i++){
	if($OP[$i]==$OP[$i-1]){
		$OP[$i]= $OP[$i];
		$strand[$i]= $strand[$i];
		$gene[$i]= "$gene[$i-1],$gene[$i]";
		$g_start[$i]=$g_start[$i-1];
		$g_end[$i]=$g_end[$i];
		$region[$i]=$region[$i];
		$eng[$i]="$eng[$i-1]:$eng[$i]";
		$eng2[$i]="$eng2[$i-1]:$eng2[$i]";
		$ratio[$i]="$ratio[$i-1]:$ratio[$i]";
#		print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\t$eng2[$i]\t$ratio[$i]\n";
	}
}

open (RatioCalc, '>Stoichiometry_ratio.csv');
my $b=0; my $wc=0;
print RatioCalc "#Operon\tStrand\tGenes\tOperon_start\tOperon_end\tdeltG_ratio\tRatio\n";
for(my $i=0; $i<$count_read1; $i++){
	if($OP[$i+1]==$OP[$i]){
		$b=1;
	}
	else{
#		print "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$region[$i]\t$eng[$i]\t$eng2[$i]\t$ratio[$i]\n";
		print RatioCalc "$OP[$i]\t$strand[$i]\t$gene[$i]\t$g_start[$i]\t$g_end[$i]\t$eng[$i]\t$ratio[$i]\n";
		$wc =$wc+1;
	}
}

close(RatioCalc);
use warnings 'uninitialized';

#my $stoi_ratio= `sed '1d' Stoichiometry_ratio.csv| wc -l`;
#print "Number of SRPS operons and their stoichiometry in the from of ratio: $stoi_ratio Stoichiometry_ratio.csv";
print "Number of SRPS operons and their stoichiometry in the from of ratio: $wc Stoichiometry_ratio.csv\n";
print "Done!\n";

system('mv tmp_slofe/rnamotif.out tmp_slofe/sequence.out tmp_slofe/rnafold.out tmp_slofe/cutoff5_stem-loops tmp_slofe/non-redundant_stem-loops tmp_slofe/final_stem-loops_sorted  .');
system('rm -rf tmp_slofe/*');
#system('mv ../rnamotif.out ../rnafold.out ../rnaM_rnaF_merged ../final_stem-loops_improved_T_intergenic ../with_eng_redef tmp_slofe/');
system('mv rnamotif.out sequence.out rnafold.out cutoff5_stem-loops non-redundant_stem-loops final_stem-loops_sorted  tmp_slofe/');

my $SSL; my $STSL;
sub system_commands{
	system('cut -f2  '.$Operon.' |sed \'1i\\\\\'|awk \'{FS=","}{print NF}\'|sed \'1d\' > tmp_slofe/operon_count');
	system('paste  '.$Operon.' tmp_slofe/operon_count > tmp_slofe/improved_operon_map1');
	system('grep -v Poly5U intergenic_intraOperonic_stem-loops |grep -v Poly4U|cut -f17|sort -g -u > tmp_slofe/tmp_srps_uniqOP');
	#Stabilizer_stem-loops-SSL
        system('grep -v Poly5U intergenic_intraOperonic_stem-loops |grep -v Poly4U|awk \'{FS="\t"}{if($19<=5)print $0}\' > Stabilizer_stem-loops-SSL');
	$SSL= `wc -l Stabilizer_stem-loops-SSL`;
	print "Number of predicted SSL: $SSL"; 
	#Stabilizer_and_terminator_stem-loops-STSL
        system('grep -v Poly5U intergenic_intraOperonic_stem-loops |grep -v Poly4U|awk \'{FS="\t"}{if($19>5)print $0}\' > Stabilizer_and_terminator_stem-loops-STSL');
	$STSL= `wc -l Stabilizer_and_terminator_stem-loops-STSL`;
	print "Number of predicted STSL: $STSL";
	system('awk \'{FS="\t"}NR==FNR{h[$1]=$2"\t"$4; next}{print $0"\t"h[$1]}\' tmp_slofe/improved_operon_map1 tmp_slofe/tmp_srps_uniqOP| cut -f2|tr \',\' \'\n\' > tmp_slofe/tmp_all_srpsGene');
	system('awk \'{FS="\t"}NR==FNR{h[$0]=1; next}{for(i in h) if($0~i){print $0"\t"i}}\' tmp_slofe/tmp_all_srpsGene tmp_slofe/improved_operon_map1|cut -f1,4,7,8 > tmp_slofe/genes_with_operon_number_new');
	system('awk \'{FS="\t"}NR==FNR{h[$3]=$1"\t"$2; next}{print $0"\t"h[$4]}\'  '.$GB.'  tmp_slofe/genes_with_operon_number_new| sort -k4|sed \'/\t1\t/d\' > tmp_slofe/srps_operon_genes_new');
	system('grep -v inside stable_stem-loops > tmp_slofe/final_stem-loops_improved_T_intergenic');
}

__END__

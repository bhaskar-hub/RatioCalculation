#!/usr/bin/perl

use strict;
use warnings;

if(!@ARGV){
	print "No FilePath Specified!\n";
	print "Usage: ./annotation_script.pl -h\n";
	exit;
}

if($ARGV[0] eq '-h' || $ARGV[0] eq '-help')
{
	help();
	exit;
}

sub help { print "This is UNIX utilities-based tool, which parses gene annotation from the genbank into a simple four-column (gene_start, gene_end, gene and strand) tab delimeted file. 
Usage: ./annotation_script.pl <annotation_file_in_genbank_format> \n Output file name will be \"parsed_gene_annotation\" \n";}

#(@ARGV == 1) || die "Usage: ./gene_annotation <annotation_file_in_genbank_format> \n Output file name will be gene_annotation_<file-name> \n";

my $op = $ARGV[0];

system_commands($op);

sub system_commands{
	#system('mkdir tmp');
	system('grep -A2 "gene    " ' .$op.'|grep -v "/db_xref\|/gene\|/pseudogene\|/old_locus_tag\|/pseudo\|CDS\|ncRNA"|sed -e \'/--/d\' -e \'s/ //g\' -e \'/^\/locus/s/$/@/\'|tr -d \'\n\'|tr \'@\' \'\n\'|sed -e \'s/gene//\' -e \'s/\.\./\t/\' -e \'s/\/locus_tag="/\t/\' -e \'s/"//\'|awk \'{FS="\t"}{if($1~/complement/) {print $1"\t"$2"\t"$3"\t-"} else {print $1"\t"$2"\t"$3"\t+"}}\'|sed -e \'s/complement(//\' -e \'s/)//\' -e \'s/>//\' -e \'s/<//\' > parsed_gene_annotation')
	}

	

	

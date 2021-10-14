# SLOFE: Stem-LOop Free-Energy
Identification of the SRPS operons and prediction of protein stoichiometry of their encoded complex using the genome sequence. SLOFE is the combination of Perl and the UNIX utilities such as awk, sed, grep etc and it runs on Linux machine.

### Table of contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Input files preparation](#Input-files-preparation)
* [USAGE](#USAGE)
* [Output files](#Output-files)
* [Error handling during RNAMotif and RNAfold installation](#Error-handling-during-RNAMotif-and-RNAfold-installation)
* [Citation](#Citation)


### Requirements

1. This tool runs on Linux machine.
2. Perl distribution v5.7.3 or above.
	> Additional Perl module: List::Utill.
3. [RNAMotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html): v3.0.7
4. [RNAfold](https://www.tbi.univie.ac.at/RNA/#download): v2.1.8  
*Note: RNAMotif and RNAfold compressed sources are provided in "dependencies" directory, as the older versions are difficult to find. It is recommended to use these sources for the installation.*

### Installation
#### Auto installation (Recommended)
1. Download the package using followng command:
	```
	git clone https://github.com/bhaskar-hub/SLOFE-RatioCalculation.git
	```
2. Installation of the RNAMotif and RNAfold from dependencies directory using the `install.sh`.
	```
	sh install.sh
	```
#### Manual installation
1. In the dependencies directory, please follow these steps to install :
	```
	tar -xvzf rnamotif-3.0.7.tar.gz
	cd rnamotif-3.0.7/
	make
	cd ../
	tar -xvzf ViennaRNA-2.1.8.tar.gz
	cd ViennaRNA-2.1.8/
	pwd
	./configure --prefix=/path/to/current/dir
	make
	```
### Input files preparation
   *Test data (test_data) and real data (Ccel_data) are given in sample_input directory. `demo.sh` file is provided in sample_input directory to run all commands.*
1. Genome of bacterial organism in fasta-format file.
2. Descriptor file for the RNAMotif. Detailed description is available in the RNAMotif [manual](http://casegroup.rutgers.edu/rnamotif.pdf). The sample descriptor file is given in sample_input directory.
3. Gene annotation file can be provided as given in sample_input directory. Alternatively, please provide gene annotation file in genbank format and run the following command:
	```
	perl scripts/annotation_script.pl gene_annotation_file_in_genbank_format
	Example:
	perl scripts/annotation_script.pl sample_input/raw_files/GCF_000008765.1_ASM876v1_genomic.gbff
	```
4. Operon map, i.e., operon annotation of the genome. SLOFE accepts operon annotation from one of these databases [Door](http://161.117.81.224/DOOR2/) or [ProOpDb](http://biocomputo2.ibt.unam.mx/OperonPredictor/). Alternatively, operon map can be provided as the same format as given in sample_input directory.
	> If operon map is from Door database, please run the following command:
	```
	perl scripts/operon_door.pl operon_map_from_Door_database
	Example:
	perl scripts/operon_door.pl sample_input/raw_files/Operon_from_door_cace
	```
	>If operon annotation is from ProOpDb, please run the following command:
	```
	perl scripts/operon_pro_op_db.pl operon_map_from_Pro_Op_Db_database
	Example:
	perl scripts/operon_pro_op_db.pl sample_input/raw_files/Operon_from_ProOpDB_cace
	```
### USAGE
#### Ratio prediction of the identified SRPS opeorns
In order to predict the ratio/stoichiometry of the SRPS operons, SLOFE uses genome_file, descriptor_file, annotation_file and operon_map. Sample input files are given in sample_input directory. Command to run SLOFE:  
*Note: `demo.sh` file is provided in sample_input directory to run all commands.*  

	perl SLOFE.pl -h
	perl SLOFE.pl <genome_fasta_file> <descr_file> <annotation_from_genbank> <operon_map> <output_directory>  
	*Note: <output_directory> parameter name should be an existing directory and it is optional. If directory is not defined, the current directory will be used for output.
	
	Example:  
	perl SLOFE.pl sample_input/test_data/test_seq.fasta sample_input/test_data/H10-mod.descr sample_input/test_data/gene_annotation sample_input/test_data/operon_map sample_input/
	perl SLOFE.pl sample_input/Ccel_data/Ccel_genome.fasta sample_input/Ccel_data/H10-mod.descr sample_input/Ccel_data/Ccel_gene_annotation sample_input/Ccel_data/Ccel_operon_map sample_input/
	
#### Scripts
1. `ratio_calc.pl` is part of `SLOFE.pl`.
2. `extract_seq.pl` extarcts the desired sequence from a fasta file using the coordinates. Usage:
	```
	perl extarct_seq.pl -h
	perl extract_seq.pl sample_input/raw_files/nc.fasta sample_input/raw_files/coords
	```
3. `maxTcontent.pl` calculates the maximum U-content and Poly(U) tails in a given sequence. Usage:
	```
	perl maxTcontent.pl -h
	perl maxTcontent.pl sample_input/raw_files/extended_right_stem_seq
	```

### Output files
**Stoichiometry_ratio.csv** : This tabulated file contains the SRPS operons and their predicted ratios. Columns are as follows:

*Operon number*

*Operon strand*

*Genes*

*Operon start*

*Operon end*

*Ratio of free-energy of stem-loops*

*Predicted Ratio*


**stable_stem-loops** : This tabulated file contains the predicted stable stem-loops from genome. Columns are as follows:

*Gene start*

*Gene end*

*Gene strand*

*Gene name*

*Stem-loop strand*

*Stem-loop start*

*Stem-loop end*

*Position of Stem-loop*

*Stem-loop length*

*Stem-loop sequence*

*Stem-loop RNA sequence*

*Stem-loop free-energy*

*Stem-loop fold*

*3'-end stem-length*

*3'-end extended 10nt sequence*

*3'-end stem plus extended 10nt sequence*

*Operon number*

*Position in operon*

*U-content of stem-loop*

*Poly(U) tail*

### Error handling during RNAMotif and RNAfold installation
1. Errors during RNAMotif installation: Possible error could be the unavailability of the "Development tools" such as "flex" and "bison". For this, either specific development tool source could be downloaded (and installed) or install using yum (`yum groupinstall "Development tools"` or `yum install flex bison`) and apt-get (`apt-get install build-essential`) based on the Linux distro.
2. Errors during RNAfold installation: `perl-devel` related errors generally occurred such as `Can't locate ExtUtils/MakeMaker.pm`. For this, either install `perl-devel` using yum/apt-get/cpan or install the specific module using yum (`yum install perl-ExtUtils-MakeMaker`) or cpan (`cpan -i ExtUtils::MakeMaker`).

### Citation
Bhaskar Y, Su X, Xu C and Xu J (2021) Predicting Selective RNA Processing and Stabilization Operons in Clostridium spp. *Front. Microbiol.* 12:673349. 
https://doi.org/10.3389/fmicb.2021.673349

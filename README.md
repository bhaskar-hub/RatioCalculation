# SLOFE: Stem-LOop Free-Energy
Identification of the SRPS operons and prediction of protein stoichiometry of their encoded complex.

### Table of contents
* [Requirements](#Requirements)
* [Installation](#Installation)
* [Input files preparation](#Input-files-preparation)
* [USAGE](#USAGE)
* [Output files](#Output-files)


### Requirements

1. This tool runs on Linux machine.
2. Perl distribution v5.7.3 or above.
	> Additional Perl module: List::Utill.
3. [RNAMotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html): v3.0.7
4. [RNAfold](https://www.tbi.univie.ac.at/RNA/#download): v2.1.8  
*Note: RNAMotif and RNAfold compressed source codes are provided in "dependencies" folder, as the older versions are difficult to find.*

### Installation
1. After installing the RNAMotif, please run `which rnamotif` if the output is `/usr/local/bin/rnamotif` or `/usr/bin/rnamotif`, then proceed to next step, otherwise, please export the path of RNAMotif folder with the following command:
	```
	export PATH=$PATH:/path/to/RNAmotif/dir
	```
	
2. After installing the RNAfold, please run `which RNAfold` if the output is `/usr/local/bin/RNAfold` or `/usr/bin/RNAfold`, then proceed to next step, otherwise, please export the path of RNAfold folder with the following command:
	```
	export PATH=$PATH:/path/to/RNAfold/dir
	```
### Input files preparation
1. Genome of bacterial organism in fasta-format file.
2. Descriptor file for the RNAMotif. Detailed description is available in the RNAMotif [manual](http://casegroup.rutgers.edu/rnamotif.pdf). The sample descriptor file is given in sample_input directory.
3. Gene annotation file can be provided as given in sample_input directory. Alternatively, please provide gene annotation file in genbank format and run the following command:
	```
	./scripts/annotation_prep gene_annotation_file_in_genbank_format
	Example:
	./scripts/annotation_prep ../sample_input/raw_files/GCF_000008765.1_ASM876v1_genomic.gbff
	```
4. Operon map, i.e., operon annotation of the genome. SLOFE accepts operon annotation from one of these databases [Door](http://161.117.81.224/DOOR2/) or [ProOpDb](http://biocomputo2.ibt.unam.mx/OperonPredictor/). Alternatively, operon map can be provided as the same format as given in sample_input directory.
	> If operon map is from Door database, please run the following command:
	```
	./scripts/operon-Door operon_map_from_Door_database
	Example:
	./scripts/operon-Door ../sample_input/raw_files/Operon_from_door_cace
	```
	>If operon annotation is from ProOpDb, please run the following command:
	```
	./scripts/operon-ProOpDb operon_map_from_Pro_Op_Db_database
	Example:
	./scripts/operon-ProOpDb ../sample_input/raw_files/Operon_from_ProOpDB_cace
	```
### USAGE
#### Ratio prediction of the identified SRPS opeorns
In order to predict the ratio/stoichiometry of the SRPS operons, SLOFE uses genome_file, descriptor_file, annotation_file and operon_map. Sample input files are given in sample_input directory. Command to run SLOFE: 

	
	./SLOFE <genome_fasta_file> <descr_file> <annotation_from_genbank> <operon_map>
	Example:
	./SLOFE sample_input/nc.fasta sample_input/H10-mod.descr sample_input/gene_annotation sample_input/operon_map
	

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

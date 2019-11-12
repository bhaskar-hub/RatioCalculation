# SLOFE: Stem-LOop Free-Energy
Identification of the SRPS operons and prediction of protein stoichiometry of their encoded complex.

### Requirements

1. This tool runs on Linux machine.
2. Perl distribution v5.7.3 or above.
	> Additional Perl module: List::Utill.
3. [RNAMotif](http://casegroup.rutgers.edu/casegr-sh-2.5.html).
4. [RNAfold](https://www.tbi.univie.ac.at/RNA/#download).

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
3. Gene annotation file as given in sample_input directory or provide gene annotation in standard genebank format.
4. Operon map, i.e., operon annotation of the genome. SLOFE accepts operon annotation from one of these databases [Door](http://161.117.81.224/DOOR2/) or [ProOpDb](http://biocomputo2.ibt.unam.mx/OperonPredictor/). Alternatively, operon map can be provided as the same format as given in sample_input directory.
	> If operon map is from Door database, please run the following command:
	```
	./scripts/operon-Door operon_map_from_Door_database
	```
	>If operon annotation is from ProOpDb, please run the following command:
	```
	./scripts/operon-Pro_Op_Db operon_map_from_Pro_Op_Db_database
	```
## Ratio prediction of the identified SRPS opeorns
In order to predict the ratio/stoichiometry of the SRPS operons, SLOFE uses genome_file, descriptor_file, annotation_file and operon_map. Sample input files are given in sample_input directory. Command to run SLOFE:

USAGE

	./SLOFE <genome_fasta-file> <descr-file> <gene_annotation> <operon_map>

## Output files
**Stoichiometry_ratio.csv** : This tabulated file contains the SRPS operons and their predicted ratios. Columns are as follows:

*Operon number*

*Operon strand*

*Genes*

*Operon start*

*Operon end*

*Ratio of free-energy of stem-loops*

*Predicted Ratio*

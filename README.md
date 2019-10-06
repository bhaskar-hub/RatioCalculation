# SLOFE: Stem-LOop Free-Energy
Identification of the SRPS operons and prediction of protein stoichiometry of their encoded complex.

### Requirements

1. This tool runs on Linux machine.
2. Perl distribution v5.7.3 or above.
	2.1 Additional perl module: List::Utill
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
3. Gene annotation file should be in standard genebank format.
4. Operon map, i.e., operon annotation of the genome should be in Door database format.
## Ratio prediction of the identified SRPS opeorns
In order to predict the ratio/stoichiometry of the SRPS operons, SLOFE uses genome_file, descriptor_file, annotation_file and operon_map. Sample input files are given in sample_input directory. Command to run SLOFE:

USAGE

	./SLOFE <genome_fasta-file> <descr-file> <gene_annotation_from_genebank> <operon_map>

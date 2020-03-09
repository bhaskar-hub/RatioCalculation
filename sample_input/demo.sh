#!bin/bash

set -e

#demo run of the SLOFE
perl ../SLOFE.pl test_data/test_seq.fasta test_data/H10-mod.descr test_data/gene_annotation test_data/operon_map

#run following command to test extarct_seq.pl
#perl ../scripts/extract_seq.pl raw_files/nc.fasta raw_files/coords > extracted_seqs

#run following command to test maxTcontent.pl
#perl ../scripts/maxTcontent.pl raw_files/extended_right_stem_seq > extended_stem_seq

#run following command to test annotation_script.pl
#perl ../scripts/annotation_script.pl raw_files/GCF_000008765.1_ASM876v1_genomic.gbff

#run following command to test operon_door.pl
#perl ../scripts/operon_door.pl raw_files/Operon_from_door_cace

#run following command to test operon_pro_op_db.pl
#perl ../scripts/operon_pro_op_db.pl raw_files/Operon_from_ProOpDB_cace

#END

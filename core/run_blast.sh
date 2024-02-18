#!/bin/bash

# dbname="df_merged_all_safe_cath_AA_seq_after_ck1-2-3.fasta"
dbname="/group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/df_merged_child_random_1000_speed_test_ForkPoolWorker-12.fasta"
##
makeblastdb -in ${dbname} -dbtype prot
(time blastp -query ${dbname} -db ${dbname} -outfmt "10 std qcovs" -out blast.csv -num_threads 32) &> blast.log

sed -i '1s/^/query,subject,%id,alignment_length,mismatches,gap_openings,query_start,query_end,subject_start,subject_end,E_value,bit_score,qcov\n/' blast.csv

rm ${dbname}.*
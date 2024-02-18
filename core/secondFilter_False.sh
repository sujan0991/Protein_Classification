#!/bin/bash

#SBATCH --job-name=scop_local_reduced # write a meaningful name
#SBATCH --time=99:00:00 # 100 hours as an example
#SBATCH --output=/group/bioinf_protstr/Ballal/TM_results_all/logs/tmv2/%x-%j.log # just write your zih name instead of sxxxxx, and keep the rest as they are 
#SBATCH --error=/group/bioinf_protstr/Ballal/TM_results_all/logs/tmv2/%x-%j.err  # just write your zih name instead of sxxxxx,  and keep the rest as they are 
#SBATCH --mem=150G # set requested memory (in MB) 
#SBATCH --ntasks=1 # leave this for now
#SBATCH --nodes=1 # leave this for n
#SBATCH --cpus-per-task=120# number of the request cpu if you have parallel processing
#SBATCH --mail-type ALL ### specify for what type of events you want to get a mail; valid options beside ALL are: BEGIN, END, FAIL, REQUEUE 


#### you must load the neede modules in advance
#### you must write module name exactly as it should be
# to see the available modules type this in terminal: module avail

## modules
module purge # to prevent conflicts 
module load apps/python3/3.9.5



python3 -c 'print("Hi everyone!")' 


echo 'Started' 

#python3 /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/Stella_project/TM_score_query_subject_v2.py 0 120
#python3 /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/core/TM_Score_calculation_v2.py 2 120
#python3 /home/bioinf/mdho200b/ballal_research_project/core/lev_distance_calculation_v3.py
#python3 /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/Stella_project/AF_against_motif_SS_local.py
python3 /group/bioinf_protstr/Ballal/Research_Project/ballal_research_project/core/local_pairwise_v3.py


echo "Finished at `date`"
echo "__DONE__"


echo "__DONE__"



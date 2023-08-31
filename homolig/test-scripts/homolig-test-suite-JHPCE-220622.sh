#Shell wrapper to run test of Homolig algorithm 
# 14 June 2022 AAG 
#$ -cwd
#$ -m e
#$ -M agirgis3@jhmi.edu
#$ -l mem_free=8G,h_vmem=8G,h_fsize=150G
#$ -pe local 1
#$ -o /users/agirgis/job-logs/
#$ -e /users/agirgis/job-logs/

module add python/3.9.10
python3 /users/agirgis/homolig/homolig/test-scripts/homolig-test-suite-JHPCE-220622.py

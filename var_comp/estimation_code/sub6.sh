#!/bin/bash
#$ -cwd
#$ -l h_data=128G,h_rt=20:00:00
#$ -o log/
#$ -e log/
#$ -N imputed
#$ -t 1-32:1
source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile

savePath=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results6/

DEBUG=true
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=2
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis6.py --tindex ${SGE_TASK_ID} --savePath $savePath

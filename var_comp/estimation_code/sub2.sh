#!/bin/bash
#$ -cwd
#$ -l h_data=96G,h_rt=10:00:00
#$ -o log/
#$ -e log/
#$ -t 1-34:1


source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile

savePath=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results2/

DEBUG=false
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=2
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis2.py --tindex ${SGE_TASK_ID} --savePath $savePath

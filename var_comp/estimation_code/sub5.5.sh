#!/bin/bash
#$ -cwd
#$ -l h_data=256G,h_rt=20:00:00
#$ -o log/
#$ -e log/
#$ -N imputed_with_VIF
#$ -t 1-33:1
source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile

savePath=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results5.5.1/

## results5.5.1: further remove the highly correlated features

DEBUG=false
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=3
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis5.5.py --tindex ${SGE_TASK_ID} --savePath $savePath

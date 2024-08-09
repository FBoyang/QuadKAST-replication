#!/bin/bash
#$ -cwd
#$ -l h_data=32G,h_rt=10:00:00
#$ -o /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/var_comp_est/realTraits/log/pc40
#$ -e /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/var_comp_est/realTraits/log/pc40
#$ -t 1-34:1
source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile

savePath=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results4/

DEBUG=false
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=32
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis4.py --tindex ${SGE_TASK_ID} --savePath $savePath

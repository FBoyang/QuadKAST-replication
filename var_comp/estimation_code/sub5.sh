#!/bin/bash
#$ -cwd
#$ -l h_data=256G,h_rt=20:00:00
#$ -o /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/var_comp_est/realTraits/log/
#$ -e /u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/var_comp_est/realTraits/log/
#$ -N imputed
#$ -t 1-1:1
source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile

savePath=/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/QuadKAST/analysis/down_stream/results5/

DEBUG=true
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=10
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis5.py --tindex ${SGE_TASK_ID} --savePath $savePath

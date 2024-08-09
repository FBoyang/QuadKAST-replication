#!/bin/bash
#$ -cwd
#$ -l h_data=32G,h_rt=10:00:00
#$ -o log/
#$ -e log/
#$ -t 1-38:1


source /u/local/Modules/default/init/modules.sh
module load gcc
module load anaconda3
source /u/home/b/boyang19/.bash_profile


# SGE_TASK_ID=10

savePath=results1.5.v2/

if [ ! -d ${savePath} ]; then
    mkdir ${savePath}
fi

DEBUG=false
if [ "$DEBUG" == true ] ; then
    SGE_TASK_ID=32
fi
echo "task ID is ${SGE_TASK_ID}"
python analysis1.5.py --tindex ${SGE_TASK_ID} --savePath $savePath

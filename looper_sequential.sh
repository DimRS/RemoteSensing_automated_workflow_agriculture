#!/bin/bash
# cd /home/stratouliasd/stars/derived/scripts/v9
# nohup ./looper_sequential.sh &>/dev/null &
export base=~/stars/derived
export baseDG=~/stars/derived/DG_v9
export basescript=~/stars/derived/scripts/v9
export R_LIBS=~/stars/rlibs

cd $baseDG/0_categ
#Rscript $basescript/0_multi_specs.R # export .csv with multi-image specifications 

# IN BASH

#DIRS=(`find . -maxdepth 1 -type d -name '*_01' | cut -c 3- `)
DIRS=( $(find . -maxdepth 1 -type d -name '*_01' | cut -c 3-) )
echo ${DIRS[142]}

for j in $(seq 0 145) # process 146 images
do
  #echo "Processing ${DIRS[j]} file...";
  $basescript/shell_command.sh ${DIRS[j]};
done
#!/bin/bash
# usage: run_project_scripts Hs inspect_search
args=("$@")
project=${args[0]}
script=${args[1]}

SOURCE=scripts/${script}.sh 
PROJ_DEST=scripts/submit/${project}_${script}.sh 
sed "s/12spe34/$project/g" $SOURCE > $PROJ_DEST 

if [ "$script" == inspect_search ] || [ "$script" == tandem_search ] 
then
	numfiles=$(ls $project/scripts/run-${script/_search/}*.sh | wc -l)
	for (( i = 0 ; i < numfiles ; i++ ))
	do
		NUM_DEST=scripts/submit/${project}_${script}${i}.sh 
		sed "s/runZnum/$i/g" $PROJ_DEST > $NUM_DEST
		qsub $NUM_DEST
	done
else
	qsub $PROJ_DEST
fi


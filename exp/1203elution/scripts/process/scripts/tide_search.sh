#!/bin/bash
#$ -V                   # Inherit the submission environment
#$ -cwd                 # Start job in submission directory
#$ -j y                 # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID
#$ -pe 4way 8   # Requests 16 tasks/node, 32 cores total
#$ -q normal
#$ -l h_rt=06:00:00     # Run time (hh:mm:ss)
#$ -M borgeson@utexas.edu
#$ -m e                # Email at Begin and End of job
#$ -P data
set -x                  # Echo commands, use "set echo" with csh
SPEC="12spe34"
SEARCH="$HOME/src.MS/tide/1.0/tide-search"
PROT="../DB/${SPEC}_all_combined_miss2.protidx"

PEP=${PROT/.protidx/.pepidx}
SUFFIX=".tideres"
cd $SPEC/tide
#$ -N tides_12spe34
for SR in $(ls *.spectrumrecords)
do
 OUT=$(basename $SR)
 OUT=${OUT/.spectrumrecords/}$SUFFIX
 $SEARCH --peptides=$PEP --proteins=$PROT --spectra=$SR --results=protobuf
 # results are stored in results.tideres
 SAVENAME="results.tideres"
 mv $SAVENAME $OUT
done

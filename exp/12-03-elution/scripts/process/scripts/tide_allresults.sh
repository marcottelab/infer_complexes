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
RESULTS="$HOME/src.MS/tide/1.0/tide-results"
PROT="$SPEC/DB/${SPEC}_all_combined_miss2.protidx"
AUX="${PROT/_miss2.protidx/.fasta.auxlocs}"
SUFFIX=".results.all"

#$ -N tider_12spe34
for SR in $(ls $SPEC/tide/*.spectrumrecords)
do
 RES=${SR/.spectrumrecords/.tideres}
 OUT=$(basename $SR)
 OUT=$SPEC/tide/${OUT/.spectrumrecords/}$SUFFIX
 $RESULTS --proteins=$PROT --spectra=$SR --results_file=$RES --out_filename=$OUT --out_format=text --aux_locations=$AUX --show_all_proteins=True
done

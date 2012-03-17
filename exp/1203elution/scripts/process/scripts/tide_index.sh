#!/bin/bash
#$ -V                   # Inherit the submission environment
#$ -cwd                 # Start job in submission directory
#$ -j y                 # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID
#$ -pe 4way 8   # Requests 16 tasks/node, 32 cores total
#$ -q largemem  # Queue name "normal"
#$ -l h_rt=08:00:00     # Run time (hh:mm:ss)
#$ -M borgeson@utexas.edu
#$ -m e                # Email at Begin and End of job
#$ -P data
set -x                  # Echo commands, use "set echo" with csh

IDX="$HOME/src.MS/crux/current/tide/tide-index"
SPEC="12spe34"
#FASTA="OMRF20110730_XENLA_EGG1_v4.mpep_trypsin_combined.fasta"
FASTA="$SPEC/DB/${SPEC}_all_combined.fasta"
PEP=${FASTA/.fasta/_miss2.pepidx}
PROT=${FASTA/.fasta/_miss2.protidx}

#$ -N tide-i-12spe34
$IDX --max_missed_cleavages=2 --enzyme=trypsin --mods_spec=C+57.02146 --fasta=$FASTA --proteins=$PROT --peptides=$PEP

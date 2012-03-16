#!/bin/bash
#$ -V                   # Inherit the submission environment
#$ -cwd                 # Start job in submission directory
#$ -j y                 # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID
#$ -pe 4way 8   # Requests 16 tasks/node, 32 cores total
#$ -q long
#$ -l h_rt=23:00:00     # Run time (hh:mm:ss)
#$ -M borgeson@utexas.edu
#$ -m e                # Email at Begin and End of job
#$ -P data
set -x                  # Echo commands, use "set echo" with csh
# replaceabales: 12spe34 for species 'Ce', runZnum for run number '2'
#$ -N insp_s_12spe34runZnum
12spe34/scripts/run-inspect-runZnum.sh

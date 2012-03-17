#!/bin/bash
#$ -V                   # Inherit the submission environment
#$ -cwd                 # Start job in submission directory
#$ -j y                 # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID
#$ -pe 4way 8
#$ -q largemem          #long, largemem
#$ -l h_rt=08:00:00     # Run time (hh:mm:ss) 8 for largemem, 24 long
#$ -M borgeson@utexas.edu
#$ -m e                # Email at Begin(b) and End(e) of job
#$ -P hpc
set -x

#$ -N tide-msconvert
/scratch/01973/blakeb/src.MS/local/bin/tide-msconvert --spectrumrecords *.mzXML

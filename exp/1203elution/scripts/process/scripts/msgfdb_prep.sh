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
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Ce/DB/Ce_all_combined.fasta -tda 0
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Dd/DB/Dd_all_combined.fasta -tda 0
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Dm/DB/Dm_all_combined.fasta -tda 0
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Hs/DB/Hs_all_combined.fasta -tda 0
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Mm/DB/Mm_all_combined.fasta -tda 0
java -Xmx10000M -cp ~/src.MS/MSGFDB/current/MSGFDB.jar msdbsearch.BuildSA -d Sp/DB/Sp_all_combined.fasta -tda 0

#! /bin/bash
abspath(){ python -c "import os.path; print os.path.abspath('$1')" ; }

usage="Usage: pepquant_setup.sh <data_root> <project_name> <output_dir>"
if [ $# -lt 3 ] ; then
    echo $usage
    exit 1
fi
# Assumes folder structure:
# data root
# - mz_ms1
#   - folders with mzXML files
# - msb_output
#   - folders with .log and _best files
# - output_dir
#   - folders with this output
# Project_name just specifies which folder.
# Must be run from inside ~/PepQuant (restriction due to finding fasta
# files--no workaround discovered yet).
# Must have a matching <species>.fasta for every 2-ch species abbrev, the first
# 2 characters of any project_name.

# constants
mzxml_dirname="mz_ms1"
msb_out_dirname="msb_output"

source_dir=$(pwd)
data_root=$1
project_name=$2
output_dir=$3
script_dir=$(dirname $(abspath $0))
pq_run=$script_dir/pepquant_ms1.py

# create project folder
project_dir=$output_dir/$project_name
if [ -d $project_dir ]; then
    echo "output/project exists:" $project_dir
    exit 1
fi
mkdir $project_dir

# link to contents
mzxml_dir=$data_root/$mzxml_dirname/$project_name
msb_out_dir=$data_root/$msb_out_dirname/$project_name
ln -s $mzxml_dir/*mzXML $project_dir
ln -s $msb_out_dir/*_best $project_dir
ln -s $msb_out_dir/*.log $project_dir

# make new sequest.params with correct fasta
fasta=${project_name:0:2}.fasta
echo "using fasta "$fasta
sed 's/FASTA_FILE/'$fasta'/g' $source_dir/sequest.params > $project_dir/sequest.params

# run pepquant
python $pq_run $fasta $project_dir

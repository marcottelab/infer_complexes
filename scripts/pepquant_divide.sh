#! /bin/bash
abspath(){ python -c "import os.path; print os.path.abspath('$1')" ; }

usage="Usage: pepquant_setup.sh <data_root> <project_name> <output_dir> <length> <run_pepquant>{true/false}"
if [ $# -lt 5 ] ; then
    echo $usage
    exit 1
fi

mzxml_dirname="mz_ms1"
msb_dirname="msb_out"
data_root=$1
project_name=$2
output_dir=$3
length=$4
do_run_pq=$5
mzxml_dir=$data_root/$mzxml_dirname/$project_name
msb_dir=$data_root/$msb_dirname/$project_name

script_dir=$(dirname $(abspath $0))
pq_setup=$script_dir/pepquant_setup.sh
pq_post=$script_dir/pepquant_post.py

n_files=$(ls $mzxml_dir/*mzXML | wc -l)
n_runs=$(expr $(expr $(expr $n_files - 1 ) / $length) + 1)
for i in $(seq $n_runs)
    do $pq_setup $data_root $project_name $output_dir $length $i $do_run_pq ;
done

if $do_run_pq; then
    project_dir=$output_dir/$project_name
    # Handles either divided or not divided output
    python $pq_post $project_dir $msb_dir $(ls -d $project_dir*)
fi

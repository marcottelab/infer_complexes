#! /bin/bash
abspath(){ python -c "import os.path; print os.path.abspath('$1')" ; }

usage="Usage: pepquant_setup.sh <data_root> <project_name> <output_dir> <length>"
if [ $# -lt 3 ] ; then
    echo $usage
    return 1
fi

mzxml_dirname="mz_ms1"
data_root=$1
project_name=$2
output_dir=$3
length=$4
mzxml_dir=$data_root/$mzxml_dirname/$project_name

script_dir=$(dirname $(abspath $0))
pq_setup=$script_dir/pepquant_setup.sh

n_files=$(ls $mzxml_dir/*mzXML | wc -l)
n_runs=$(expr $(expr $n_files / $length) + 1)
for $index in $(seq $n_runs)
    do $pq_setup $data_root $project_name $output_dir $length $index &
done
wait

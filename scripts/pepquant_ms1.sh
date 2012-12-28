#! /bin/bash
usage="Usage: pepquant_ms1.sh <fasta_in_source_dir> <project_dir>
<python_ms1_script>"
args=("$@")
if [ ${#args[@]} -lt 3 ]; then
    echo $usage
    exit 1
fi

fasta=$1
source_dir=$(dirname $fasta)
project_dir=$2
py_script=$3

java -jar $source_dir/MSGFplusExtractor.jar $project_dir
python $py_script $fasta $source_dir/seq_header.txt $project_dir/*_best
$source_dir/genPepList.pl $project_dir
$source_dir/runIPCbatch.pl $project_dir/pep4quant_sequest.txt
$source_dir/scanDir4Quant.pl $project_dir $project_dir/TOTAL.ISO 0.5 10
for f in $(ls -d $project_dir/*mzXML_dta); do {
    cmd="java -Xmx3000M -cp $source_dir/PepQuant.jar scanDir $f
    $project_dir/ISO_QUAN $project_dir/TOTAL.ISO 0.5 10"
    echo $cmd;
    $cmd;
} done
$source_dir/genMS1MS2RTmap.pl $project_dir
$source_dir/pepQuan.pl $project_dir/ISO_QUAN 0.99

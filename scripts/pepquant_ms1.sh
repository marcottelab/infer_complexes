#! /bin/bash
abspath(){ python -c "import os.path; print os.path.abspath('$1')" ; }

usage="Usage: pepquant_ms1.sh <fasta_in_source_dir> <project_dir>"
if [ $# -lt 2 ] ; then
    echo $usage
    exit 1
fi
fasta=$1
source_dir=$(dirname $fasta)
project_dir=$2
py_dir=$(dirname $(abspath $0))
py_convert=$py_dir/msgf2seq.py
py_combine=$py_dir/sequest2txt.py

# test python version
echo "python:" $(which python)
# extract mzXML files to mzXML_data directories
java -Xmx500M -jar $source_dir/MSGFplusExtractor.jar $project_dir
# convert output from search engines and msblender filter into sequest format
python $py_convert $fasta $project_dir/*pep_count_FDR*.log $project_dir/*_best
# combine output files from separate search engines into sequest input
python $py_combine $source_dir/seq_header.txt $project_dir/*.sequestformat
# run PepQuant
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

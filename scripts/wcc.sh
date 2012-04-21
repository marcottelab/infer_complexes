#! /bin/bash
usage="Usage: wcc.sh <filename> <width>"
# run from where you want target directory made
args=("$@")
if [ ${#args[@]} -lt 2 ]; then
    echo $usage
    exit 1
fi

fname=${args[0]}
width=${args[1]}
d=$(dirname $0)
fname_out=$fname.T

#python -c "import sys; for l in sys.stdin.readlines(): if l.strip(): print '\t'.join(zip(*(l.split('\t'))))" < $fname > $fname_out
#python -c "import sys;  print ['\t'.join(zip(*(l.split()))) for l in sys.stdin.readlines() if l.strip() and l[0]!='#']" < $fname > $fname_out
python $d/transpose.py $fname $fname_out


Rscript $d/_wcc.R $fname_out $width 
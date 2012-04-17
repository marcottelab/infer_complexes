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

python -c "import sys; print('\n'.join('\t'.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $fname > $fname_out

Rscript $d/_wcc.R $fname_out $width 
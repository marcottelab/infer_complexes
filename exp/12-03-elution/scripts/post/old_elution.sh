#! /bin/bash
# usage: run.sh sourcedir match 
# run from within destination directory
msb="$HOME/Dropbox/complex/tools/MSblender"
fdr_string=005
fdr_num=0.05
args=("$@")
sourcedir=${args[0]}
match=${args[1]}

for bestfile in $(ls $sourcedir/* | grep ${match})
do
	filebase=$(basename "$bestfile")
	filebase="${filebase%%.*}"
	echo "filebase: $filebase"
	sed "s/^000/${filebase}/g" $sourcedir/$bestfile > $filebase.best
done

echo "combining into $match.best"
# put the header lines from the first file in a combined .best file
head -n 2 $(ls *.best | grep ${match} | head -n 1) > $match.combined_best
# put all contents except the first line from all files into combined file
for awkfile in $(ls *.best | grep ${match})
do
    awk 'FNR>2' $awkfile >> $match.combined_best
done 
# -e enables special characters
echo -e "Tide\t${match}.combined_best" >> $match.conf
echo "creating msblender_in"
$msb/pre/make-msblender_in.py $match.conf 
echo "running msblender"
$msb/src/msblender $match.msblender_in 
$msb/post/make-spcount.py $match.msblender_in.msblender_out $match.prot_list $fdr_num
$msb/post/filter-msblender.py $match.msblender_in.msblender_out $fdr_string > $match.filter
tail -n 1 $match.filter

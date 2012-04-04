#! /bin/bash
# usage: elution.sh sourcedir shortname match 
# run from where you want target directory made
#msb="$HOME/git/MSblender" # on tacc
msb="$HOME/Dropbox/complex/tools/MSblender" # on lab network
fdr_string=001
fdr_num=0.01
args=("$@")
sourcedir=${args[0]}
shortname=${args[1]}
match=${args[2]}
searches=( Tide Inspect MSGFDB )
extensions=( xcorr_hit_list_best MQscore_hit_list_best logSpecProb_hit_list_best )

if [ -d $shortname ]; then
    echo "directory exists."
    exit 1
fi

mkdir $shortname
cd $shortname
echo "sourcedir:"$sourcedir", combining into $shortname.(search).combined_best"
nsearches=${#searches[@]}
for (( i=0;i<$nsearches;i++ ))
do
	search=${searches[$i]}
	ext=${extensions[$i]}
	echo "match"${match}*.$ext
	comb_file=$shortname.$search.combined_best
	for bestfile in $(ls $sourcedir/${match}*.$ext)
	do
		if [ $search != Tide ]; then
			cat $bestfile >> $shortname.$search.combined_best
		else
			# Tide is missing spectrum ids--handle separately
			filebase=$(basename "$bestfile")
			filebase="${filebase%%.*}"
			echo "filebase: $filebase"
			sed "s/^000/${filebase}/g" $bestfile >> $comb_file
		fi
	done
	# -e enables special characters, needed for tab and maybe newline too
	echo -e "${search}\t${comb_file}" >> $shortname.conf
done

echo "creating msblender_in"
$msb/pre/make-msblender_in.py $shortname.conf 
echo "running msblender"
$msb/src/msblender $shortname.msblender_in 
$msb/post/make-spcount.py $shortname.msblender_in.msblender_out $shortname.prot_list $fdr_num
$msb/post/filter-msblender.py $shortname.msblender_in.msblender_out $fdr_string > $shortname.filter
tail -n 1 $shortname.filter

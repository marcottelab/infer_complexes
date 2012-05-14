#! /bin/bash
usage="Usage: weka.sh train_key elut_file_pattern data_scores species fnet_gene_dict suffix "
flog='logtest.txt'
args=("$@")
if [ ${#args[@]} -lt 5 ]; then
    echo $usage
    exit 1
fi

train_key=$1
species=$4
suffix=$6
d=$(dirname $0)
python $d/../ml.py $train_key $2 $3 $species $5 $suffix
fname=${species}_${train_key}_$suffix
ftrain=${fname}_train.arff
ftest=${fname}_test.arff
fwekaout=${fname}_scored.tab
ffilt=${fname}_filtered.tab
#-cp /Users/blakeweb/Dropbox/complex/tools/weka-3-6-6/weka.jar 
java weka.classifiers.meta.LogitBoost -P 100 -F 0 -R 1 -L -1.7976931348623157E308 -H 1.0 -S 1 -I 10 -W weka.classifiers.trees.DecisionStump -t $ftrain -T $ftest -p 0 > ${fwekaout}

# First 5 lines of weka output are header
echo $(wc -l $fwekaout)" output predictions from weka"
./$d/weka_filter.sh $fwekaout $ffilt
python $d/../cv.py $ffilt | tee -a $flog
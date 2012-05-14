#! /bin/bash
fscored=$1
ffilt=$2
# Look for 'true' in the third column, which means it will either follow true, false, or ? (if it's not a test set).  Get rid of '1:' and '2:', class labels.  Get rid of the error column which contains a '+' or nothing. Get just columns 4, 2, and 3.  Replace space separation with tab separation.  Sort descending by the score, which is now first.
sed '1,5d' $fscored | grep -e"[e\?][ ]*1:true" | sed 's/1\://g' | sed 's/2\://g' | sed 's/+//g' | awk '{print $4,$2,$3}' | tr ' ' '\t' | sort -r > $ffilt
echo $(wc -l $ffilt)" interactions with posterior>0.5 in "${ffilt}
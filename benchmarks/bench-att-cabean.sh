#!/bin/bash
export LC_NUMERIC="en_US.UTF-8"

CURR=`pwd`
MODEL_DIR = $1
TIME_LIMIT = $2
INPUTDIR=$CURR/$MODEL_DIR

rm -f $INPUTDIR/*.ispl
for j in `ls $INPUTDIR`
do
	echo "Run $j"
	#start=$(date +%s.%N)
	second="s"
	timeout_str="$timeout$second"
	timeout $timeout_str ./cabean -compositional 2 $INPUTDIR/$j;killall -9 cabean
	#dur=$(echo "$(date +%s.%N) - $start" | bc)
	#execution_time=`printf "%.2f seconds" $dur`
	#sed -i '1s/^/'"$execution_time"'\n/' $OUTDIR/$j
	echo "Finish $j"
done

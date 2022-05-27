#!/bin/bash

BASEDIR=.ksd
METHOD=AIC

for dir_no in `seq 1 5`;
do
	dir="Melissa_Model${dir_no}"
        for rep in `seq 0 9`;
        do
		file="model${dir_no}_replicate${rep}.stru"
		for k in 1 2 3 4
		do
			VAL=`grep $METHOD $BASEDIR/$dir/${file}.admix.K=${k}.out.txt | cut -d ' ' -f 3`
			if [ $k == 1 ]; then
				MIN=$VAL
				MINK=$k
			elif (( $(echo "$VAL < $MIN" | bc -l) )); then
				MIN=$VAL
				MINK=$k
			fi
		done
		echo "$file $rep $MINK"
	done
done

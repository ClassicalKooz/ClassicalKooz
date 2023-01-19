#!/bin/sh

for i in 1 2 3 4 5 6
do
	cd Melissa_Model$i
	for j in 0 1 2 3 4 5 6 7 8 9
	do
		sed -i '1d' replicate$j.str.recode.strct_in
		mv replicate$j.str.recode.strct_in replicate$j.str
	done
	cd ..
done


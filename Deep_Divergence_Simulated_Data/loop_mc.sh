#!/bin/sh

for i in 1 2 3 4 5 6
do
	cd Melissa_Model$i
	for j in 0 1 2 3 4 5 6 7 8 9
	do
		sed -i 's/tsk//g' replicate$j.str
		for k in 1 2 3 4
		do
			multiclust -f replicate$j.str -a -k $k -s 4
		done
	done
	cd ..
done

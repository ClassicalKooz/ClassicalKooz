#!/bin/sh

for i in 1 2 3 4 5
do
	cd Model$i
	for j in 1 2 3 4 5 6 7 8 9 10
	do
		sed -i 's/tsk//g' model$i_replicate$j.str
		for k in 1 2 3 4
		do
			multiclust -f model$i_replicate$j.str -a -k $k -s 4
		done
	done
	cd ..
done

#!/bin/bash


for j in {1..10}
do
	cat header.txt >> rep"$j".u
	for i in {1..20}
	do
		locuslength=$(head -1 model2_replicate"$j"_"$i".u | awk '{print $2}' | wc -c)
		echo locus"$i" 10 10 10 $locuslength I 1  >> rep"$j".u
		cat model2_replicate"$j"_"$i".u >> rep"$j".u
	done
done

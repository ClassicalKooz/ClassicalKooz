#!/bin/sh

for i in 0 1 2 3 4 5 6 7 8 9
	do
		ls model1_replicate$i*.vcf > replicate$i
		bcftools concat -f replicate$i -o replicate$i.vcf -O v
	done



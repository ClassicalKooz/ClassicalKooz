#!/bin/sh

for i in 0 1 2 3 4 5 6 7 8 9
do
	sed -i 's/tsk_/tsk/g' replicate$i.vcf
	plink --vcf replicate$i.vcf --recode structure --out replicate$i.str
done

#!/usr/bin/env bash

for FILE in *.vcf;
do echo $FILE;
vcftools --vcf $FILE --weir-fst-pop 'pop1.txt' --weir-fst-pop 'pop2.txt' --weir-fst-pop 'ghost.txt' --fst-window-size 100 --fst-window-step 100 --out $FILE.3pop.fst
done
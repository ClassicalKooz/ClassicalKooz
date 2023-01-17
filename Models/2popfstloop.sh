#!/usr/bin/env bash
for file in /Models;
do
for i in {1..15};
    do 
    echo Model$i
    for *.vcf in Model$1;
        do echo $FILE;
        done
#vcftools --vcf $FILE --weir-fst-pop 'pop1.txt' --weir-fst-pop 'pop2.txt' --fst-window-size 100 --fst-window-step 100 --out $FILE.2pop.fst
    done
done

for i in {1..10};
do cd replicate$i
    for FILE in *.vcf;
    do echo $FILE;
    vcftools --vcf $FILE --weir-fst-pop pop1 --weir-fst-pop pop2 --weir-fst-pop ghost --fst-window-size 100 --fst-window-step 100 --out $FILE.3pop.fst
    done
cd ..
done
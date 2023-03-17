for i in {1..10};
do cd replicate$i
    for FILE in *.vcf;
    do echo $FILE;
    vcftools --vcf $FILE --window-pi 100 --window-pi-step 100 --keep pop1 --keep pop2 --out $FILE.2pop.pi
    done
cd ..
done
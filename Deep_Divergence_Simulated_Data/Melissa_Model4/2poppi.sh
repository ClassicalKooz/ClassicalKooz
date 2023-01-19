for FILE in *.vcf;
do echo $FILE;
vcftools --vcf $FILE --window-pi 100 --window-pi-step 100 --keep pop1.txt --keep pop2.txt --out $FILE.2pop.pi
done
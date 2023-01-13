for FILE in *.vcf;
do echo $FILE;
vcftools --vcf $FILE --window-pi 100 --window-pi-step 100 --keep pop1.txt --keep pop2.txt --keep ghost.txt --out $FILE.3pop.pi
done
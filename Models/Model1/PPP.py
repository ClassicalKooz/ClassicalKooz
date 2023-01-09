import sys
import os
import subprocess

from pgpipe import vcf_split_pysam, vcf_to_ima, vcf_filter, vcf_calc, vcf_phase, stat_sampler, vcf_split
#took out four_gamete and sampler
from pgpipe.logging_module import initLogger
from pgpipe.informative_loci_filter import filter_bed_regions
#from pgpipe.subtract_bed import filter_stat
#causes error "no module named pgpipe.substract_bed
import pysam
import pgpipe.vcf_filter as vcf_filter

print ("Imports complete")

rootdir = '/home/faculty/asethuraman/Documents/MichaelK/Models'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if 'vcf' in file:
            print(os.path.join(subdir, file))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="TajimaD",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="TajimaD",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt", model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
            vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="/home/faculty/asethuraman/Documents/MichaelK/models.model",model="2Pop",out_prefix="2pop")





#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="TajimaD",model_file="models.txt",model="2Pop",statistic_window_size=100)
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="TajimaD",model_file="models.txt",model="3Pop",statistic_window_size=100)
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="het-fit",model_file="models.txt",model="2Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="het-fit",model_file="models.txt",model="3Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="2Pop",out_prefix="2pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="3Pop",out_prefix="3pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="2Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="3Pop", out_prefix="3pop")



 

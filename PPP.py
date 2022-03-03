
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
import vcf_to_sfs as vs
#import vcf_format_conversions as vfc

print ("Imports complete")

rootdir = '/home/faculty/asethuraman/Documents/MichaelK/Models'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if 'vcf' in file:
            if os.path.exists("/home/faculty/asethuraman/Documents/MichaelK/Models/Statistic_Files"):
                os.rmdir("/home/faculty/asethuraman/Documents/MichaelK/Models/Statistic_Files")
                
                
           
            print(os.path.join(subdir, file))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="TajimaD",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="TajimaD",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt", model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="/home/faculty/asethuraman/Documents/MichaelK/models.model",model="2Pop",out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="/home/faculty/asethuraman/Documents/MichaelK/models.model",model="3Pop",out_prefix=os.path.join(subdir,"3pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="/home/faculty/asethuraman/Documents/MichaelK/models.model",model="2Pop",out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="/home/faculty/asethuraman/Documents/MichaelK/models.model",model="3Pop",out_prefix=os.path.join(subdir,"3pop"))
            #vs.build_sfs(os.path.join(subdir, file),"/home/faculty/asethuraman/Documents/MichaelK/models.model","2Pop",downsamplesizes=[10,10],folded="FALSE",out=os.path.join(subdir,"2pop.sfs"))
            #vs.build_sfs(os.path.join(subdir, file),"/home/faculty/asethuraman/Documents/MichaelK/models.model","3Pop",downsamplesizes=[10,10,10],folded="FALSE",out=os.path.join(subdir,"3pop.sfs"))
            #(vcf= os.path.join(subdir, file), out_format="ped",out_prefix=os.path.join(subdir,"PED"))
            #vcf_format_conversions(vcf= "Model1/replicate1/model1_replicate1.vcf", out_format="ped-12",out_prefix=os.path.join(subdir,"PED"))
            



#programming note: make an thing that removes "statistic files" before moving on to the next window-pi, should solve it
#follow up^: that worked and window-pi 3 pop is good to go
#02/16/22 all tests completed
#03/02/22 all test completed but with binary           

#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="TajimaD",model_file="models.txt",model="2Pop",statistic_window_size=100)
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="TajimaD",model_file="models.txt",model="3Pop",statistic_window_size=100)
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="het-fit",model_file="models.txt",model="2Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="het-fit",model_file="models.txt",model="3Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="2Pop",out_prefix="2pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="3Pop",out_prefix="3pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="2Pop")
#vcf_calc.run(vcf="/home/faculty/asethuraman/Documents/MichaelK/example.vcf",calc_statistic="window-pi",statistic_window_size=100,statistic_window_step=100,model_file="models.model",model="3Pop", out_prefix="3pop")

print ("Done!")
 

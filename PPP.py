import sys
import os
import subprocess

from pgpipe import four_gamete, vcf_split_pysam, vcf_to_ima, vcf_filter, vcf_calc, vcf_sampler, vcf_phase, stat_sampler, vcf_split
from pgpipe.logging_module import initLogger
from pgpipe.informative_loci_filter import filter_bed_regions
from pgpipe.subtract_bed import filter_stat
import pysam

print ("Imports complete")

path="models2/model1"
dir_list = os.listdir(path)
#print(dir_list)
for root, dirs, files in os.walk(path, topdown=False):
   for name in files:
        #if 'vcf' in name:
            #print(os.path.join(root, name))
          
        #if os.path.exists("model2/model1"):
                #os.rmdir("Melissa_Simulated_Data/Statistic_Files") 
            #vcf_calc.run(vcf = os.path.join(root, name),calc_statistic="TajimaD",model_file="Models",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(root,"2pop"))
            #vcf_calc.run(vcf=os.path.join(subdir, file),calc_statistic="TajimaD",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt", model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="2Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"2pop"))
            #vcf_calc.run(vcf = os.path.join(subdir, file),calc_statistic="het-fit",model_file="/home/faculty/asethuraman/Documents/MichaelK/models.txt",model="3Pop",statistic_window_size=100,out_prefix=os.path.join(subdir,"3pop"))
<<<<<<< HEAD
        vcf_calc.run(vcf=os.path.join(root, name),calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="/home/michael/GhostPopulationSimulation/models.model",model="2Pop",out_prefix=os.path.join(name,"2pop.fst"))
=======
            #vcf_calc.run(vcf=os.path.join(root, name),calc_statistic="windowed-weir-fst",statistic_window_size=100,statistic_window_step=100,model_file="/home/michael/GhostPopulationSimulation/models.model",model="2Pop",out_prefix=os.path.join(name,"2pop.fst"))
>>>>>>> 227be4475bc73163c7778c0cc36b007033c19add
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
 
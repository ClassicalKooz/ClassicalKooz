from ast import Num
from pickle import FALSE, TRUE
import msprime
import stdpopsim
import random
import sys
import numpy as np
from IPython.display import SVG
import pandas as pd

#Specifications of the demographic model
demography = msprime.Demography.island_model([100]*3,migration_rate=0.001)

#TODO for MK: We need to write several loops. 1) We have to loop over migration rates (0 - no gene flow, 0.1 - low gene flow, 0.5 - high gene flow)
# Keep divergence times constant along the lines of divergences of modern humans from Neanderthals - (1850 generations - recent, 18500 generations - ancient)
#Keep population sizes constant - keep all pops at Ne = 5000
#save the graphs and the VCF files
#do at least 10 replicates/model, each replicate with a different random number seed
Fst34=[]
Fst32=[]
Fst42=[]
Divergence34=[]
Divergence32=[]
Divergence42=[]
Diversity3=[]
Diversity4=[]
Diversity2=[]
TajimaD=[]
#write a for loop over the next three commands (1) simulate tree, (2) simulate mutations, (3) write a vcf
for k in range(0,10): # looping over 10 replicates
    for j in range(0,20): # looping over number of independent loci (20)
        x=random.randint(1,9999)
        #Call msprime to simulate tree sequence under the demographic model
        ts = msprime.sim_ancestry(sequence_length=1000, samples={"pop_0": 10, "pop_1": 10,"pop_2":10,"pop_3":10,"pop_4":10}, demography=demography, random_seed=x)
        #Simulate mutations along the tree sequence simulated by msprime above
        ts = msprime.sim_mutations(ts, rate=0.000001, random_seed=x)#, model='binary')

        #Write a vcf file containing all the mutations simulated under the model above - example.vcf should now contain the variants
        with open("ILM3_replicate%r_%r.vcf" %(k,j), "w") as vcf_file:
            ts.write_vcf(vcf_file,contig_id=str(j))
            ts.write_fasta("ILM3_replicate%r_%r.fasta" %(k,j))
            imfile = open("ILM3_replicate%r_%r.u" %(k,j),"w",buffering=1)
            for i,h in enumerate(ts.haplotypes()):
                print(f"Sample{i} {h}",file=imfile)
        #calculate all stats
        
        Fst34.append(ts.Fst(sample_sets=[ts.samples(population=0), ts.samples(population=1)],mode="site"))
        Fst32.append(ts.Fst(sample_sets=[ts.samples(population=0), ts.samples(population=2)],mode="site"))
        Fst42.append(ts.Fst(sample_sets=[ts.samples(population=1), ts.samples(population=2)],mode="site"))
        
        Divergence34.append(ts.divergence(sample_sets=[ts.samples(population=0), ts.samples(population=1)],mode="site"))
        Divergence32.append(ts.divergence(sample_sets=[ts.samples(population=0), ts.samples(population=2)],mode="site"))
        Divergence42.append(ts.divergence(sample_sets=[ts.samples(population=1), ts.samples(population=2)],mode="site"))

        Diversity3.append(ts.diversity(sample_sets=ts.samples(population=0),mode="site"))
        Diversity4.append(ts.diversity(sample_sets=ts.samples(population=1),mode="site"))
        Diversity2.append(ts.diversity(sample_sets=ts.samples(population=2),mode="site"))
        
        TajimaD.append(ts.Tajimas_D(sample_sets=[ts.samples(population=0), ts.samples(population=1),ts.samples(population=2)],mode="site"))
        
        
        
'''d={"Fst34":Fst34,"Fst32":Fst32,"Fst42":Fst42}
df=pd.DataFrame(data=d)
df.to_csv('SS3_Fst.csv',index=False)
d={"Divergence34":Divergence34,"Divergence32":Divergence32,"Divergence42":Divergence42}
df=pd.DataFrame(data=d)
df.to_csv('SS3_Divergence.csv',index=False)
Diversity=[Diversity3,Diversity4,Diversity2]
df=pd.DataFrame(data=Diversity)
df=df.transpose()
df.to_csv('SS3_Diversity.csv',index=False)
df=pd.DataFrame(data=TajimaD)
df.to_csv('SS3_TajimaD.csv',index=False)'''

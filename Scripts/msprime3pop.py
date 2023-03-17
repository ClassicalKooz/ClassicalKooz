from ast import Num
from pickle import FALSE, TRUE
import msprime
import random
import sys
import numpy as np
from IPython.display import SVG
import pandas as pd

#Specifications of the demographic model
demography = msprime.Demography()
demography.add_population(name="NA1", initial_size=5_000)
demography.add_population(name="NA2", initial_size=5_000)
# Ghost is ghost
demography.add_population(name="Ghost", initial_size=5_000)
demography.add_population(name="pop1", initial_size=5_000)
demography.add_population(name="pop2", initial_size=5_000)
demography.add_population_split(time=18500, derived=["NA2", "Ghost"], ancestral="NA1")
demography.add_population_split(time=1850, derived=["pop1", "pop2"], ancestral="NA2")
demography.set_migration_rate(source="pop1", dest="pop2", rate=0.010)
demography.set_migration_rate(source="pop2", dest="pop1", rate=0.010)
demography.set_migration_rate(source="pop1", dest="Ghost", rate=0.000)
demography.set_migration_rate(source="Ghost", dest="pop1", rate=0.050)
demography.set_migration_rate(source="pop2", dest="Ghost", rate=0.000)
demography.set_migration_rate(source="Ghost", dest="pop2", rate=0.050)
demography.sort_events()

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
for k in range(1,11): # looping over 10 replicates
    for j in range(1,21): # looping over number of independent loci (20)
        x=random.randint(1,9999)
        #Call msprime to simulate tree sequence under the demographic model
        ts = msprime.sim_ancestry(sequence_length=1000, samples={"pop1": 10, "pop2": 10, "Ghost":10}, demography=demography, random_seed=x)
        #Simulate mutations along the tree sequence simulated by msprime above
        #x=random.randint(1,9999)
        ts = msprime.sim_mutations(ts, rate=0.00000001, random_seed=x)#, model='binary')

        #Write a vcf file containing all the mutations simulated under the model above - example.vcf should now contain the variants
        with open("model2_replicate%r_%r.vcf" %(k,j), "w") as vcf_file:
            #ts.write_vcf(vcf_file,contig_id=str(j))
            ts.write_fasta("model5_replicate%r_%r.fasta" %(k,j))
            imfile = open("model5_replicate%r_%r.u" %(k,j),"w",buffering=1)
            for i,h in enumerate(ts.haplotypes()):
                print(f"Sample{i} {h}",file=imfile)
        #calculate all stats
        
        #Fst34.append(ts.Fst(sample_sets=[ts.samples(population=3), ts.samples(population=4)],mode="site"))
        #Fst32.append(ts.Fst(sample_sets=[ts.samples(population=3), ts.samples(population=2)],mode="site"))
        #Fst42.append(ts.Fst(sample_sets=[ts.samples(population=4), ts.samples(population=2)],mode="site"))
        '''
        Divergence34.append(ts.divergence(sample_sets=[ts.samples(population=3), ts.samples(population=4)],mode="site"))
        Divergence32.append(ts.divergence(sample_sets=[ts.samples(population=3), ts.samples(population=2)],mode="site"))
        Divergence42.append(ts.divergence(sample_sets=[ts.samples(population=4), ts.samples(population=2)],mode="site"))

        Diversity3.append(ts.diversity(sample_sets=ts.samples(population=3),mode="site"))
        Diversity4.append(ts.diversity(sample_sets=ts.samples(population=4),mode="site"))
        Diversity2.append(ts.diversity(sample_sets=ts.samples(population=2),mode="site"))
        
        TajimaD.append(ts.Tajimas_D(sample_sets=[ts.samples(population=3), ts.samples(population=4),ts.samples(population=2)],mode="site"))
        '''
        
        
#d={"Fst34":Fst34,"Fst32":Fst32,"Fst42":Fst42}
#df=pd.DataFrame(data=d)
#df.to_csv('Model1_Fst.csv',index=False)
'''d={"Divergence34":Divergence34,"Divergence32":Divergence32,"Divergence42":Divergence42}
df=pd.DataFrame(data=d)
df.to_csv('Model5_Divergence.csv',index=False)
Diversity=[Diversity3,Diversity4,Diversity2]
df=pd.DataFrame(data=Diversity)
df=df.transpose()
df.to_csv('Model5_Diversity.csv',index=False)
df=pd.DataFrame(data=TajimaD)
df.to_csv('Model6_TajimaD.csv',index=False)'''




#Prints the demography 
#print(demography)

#Graphical output of the simulated demographic model
#import demesdraw
#import matplotlib.pyplot as plt 
#graph = msprime.Demography.to_demes(demography)
#fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
#demesdraw.tubes(graph, ax=ax, seed=1)
#plt.show() 




from pickle import TRUE
import msprime
import stdpopsim
import random
import sys

#Specifications of the demographic model
demography = msprime.Demography()
demography.add_population(name="NA1", initial_size=2_500)
demography.add_population(name="NA2", initial_size=2_500)
# Ghost is ghost
demography.add_population(name="Ghost", initial_size=2_500)
demography.add_population(name="pop1", initial_size=2_500)
demography.add_population(name="pop2", initial_size=2_500)
demography.add_population_split(time=2000000, derived=["NA2", "Ghost"], ancestral="NA1")
demography.add_population_split(time=500000, derived=["pop1", "pop2"], ancestral="NA2")
demography.set_migration_rate(source="pop1", dest="pop2", rate=0.000)
demography.set_migration_rate(source="pop2", dest="pop1", rate=0.000)
demography.set_migration_rate(source="pop1", dest="Ghost", rate=0.00)
demography.set_migration_rate(source="Ghost", dest="pop1", rate=0.00)
demography.set_migration_rate(source="pop2", dest="Ghost", rate=0.00)
demography.set_migration_rate(source="Ghost", dest="pop2", rate=0.00)
demography.sort_events()

#TODO for MK: We need to write several loops. 1) We have to loop over migration rates (0 - no gene flow, 0.1 - low gene flow, 0.5 - high gene flow)
# Keep divergence times constant along the lines of divergences of modern humans from Neanderthals - (1850 generations - recent, 18500 generations - ancient)
#Keep population sizes constant - keep all pops at Ne = 5000
#save the graphs and the VCF files
#do at least 10 replicates/model, each replicate with a different random number seed

#write a for loop over the next three commands (1) simulate tree, (2) simulate mutations, (3) write a vcf
for k in range(0,1): # looping over replicates
    for j in range(0,10): # looping over number of independent loci
        x=random.randint(1,9999)
        #Call msprime to simulate tree sequence under the demographic model
        ts = msprime.sim_ancestry(sequence_length=1000, samples={"pop1": 10, "pop2": 10, "Ghost":10}, demography=demography, random_seed=x)
        #Simulate mutations along the tree sequence simulated by msprime above
        mts = msprime.sim_mutations(ts, rate=0.000001, random_seed=x)#, model='binary')

        #Write a vcf file containing all the mutations simulated under the model above - example.vcf should now contain the variants
        with open("model6_replicate%r_%r.vcf" %(k,j), "w") as vcf_file:
            mts.write_vcf(vcf_file,contig_id=str(j))
        #mts.write_fasta("model6_replicate%r_%r.fasta" %(i,j))
            imfile = open("model6_replicate%r_%r.u" %(k,j),"w",buffering=1)
            for i,h in enumerate(mts.haplotypes()):
                print(f"Sample{i} {h}",file=imfile)




#Prints the demography 
#print(demography)

#Graphical output of the simulated demographic model
#import demesdraw
#import matplotlib.pyplot as plt 
#graph = msprime.Demography.to_demes(demography)
#fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
#demesdraw.tubes(graph, ax=ax, seed=1)
#plt.show() 




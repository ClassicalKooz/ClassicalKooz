import msprime
import stdpopsim

#Specifications of the demographic model
demography = msprime.Demography()
demography.add_population(name="NA1", initial_size=5_000)
demography.add_population(name="NA2", initial_size=5_000)
# Ghost is ghost
demography.add_population(name="Ghost", initial_size=5_000)
demography.add_population(name="pop1", initial_size=5_000)
demography.add_population(name="pop2", initial_size=5_000)
demography.add_population_split(time=18500, derived=["pop2", "NA2"], ancestral="NA1")
demography.add_population_split(time=1850, derived=["Ghost", "pop1"], ancestral="NA2")

demography.set_migration_rate(source="pop2", dest="Ghost", rate=0.1)
demography.set_migration_rate(source="Ghost", dest="pop2", rate=0)
demography.set_migration_rate(source="pop1", dest="pop2", rate=0.1)
demography.set_migration_rate(source="pop2", dest="pop1", rate=0.1)
demography.set_migration_rate(source="pop1", dest="Ghost", rate=0.1)
demography.set_migration_rate(source="Ghost", dest="pop1", rate=0)

demography.sort_events()

#TODO for MK: We need to write several loops. 1) We have to loop over migration rates (0 - no gene flow, 0.1 - low gene flow, 0.5 - high gene flow)
# Keep divergence times constant along the lines of divergences of modern humans from Neanderthals - (1850 generations - recent, 18500 generations - ancient)
#Keep population sizes constant - keep all pops at Ne = 5000
#save the graphs and the VCF files
#do at least 10 replicates/model, each replicate with a different random number seed

#write a for loop over the next three commands (1) simulate tree, (2) simulate mutations, (3) write a vcf

#Call msprime to simulate tree sequence under the demographic model
ts = msprime.sim_ancestry(sequence_length=100, samples={"pop1": 10, "pop2": 10, "Ghost":10}, demography=demography, random_seed=9695)

#Simulate mutations along the tree sequence simulated by msprime above
mts = msprime.sim_mutations(ts, rate=0.01, random_seed=9695)

#Write a vcf file containing all the mutations simulated under the model above - example.vcf should now contain the variants
with open("model11_replicate3.vcf","w") as vcf_file:
    mts.write_vcf(vcf_file,contig_id="0")

#Prints the demography 
print(demography)

#Graphical output of the simulated demographic model
import demesdraw
import matplotlib.pyplot as plt 
graph = msprime.Demography.to_demes(demography)
fig, ax = plt.subplots()  # use plt.rcParams["figure.figsize"]
demesdraw.tubes(graph, ax=ax, seed=1)
#plt.show() 




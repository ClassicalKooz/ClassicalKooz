import seaborn
import pandas as pd
df=pd.read_csv('/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/Model1_Fst_3pop_combined')
violin = seaborn.violinplot(data=df, x="MEAN_FST")
violin
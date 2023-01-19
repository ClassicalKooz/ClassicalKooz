import seaborn
import pandas as pd
df=pd.read_csv('/home/michael/GhostPopulationSimulation/Models/All_Models_Fst_combined')
violin = seaborn.violinplot(data=df, x="Model",y="MEAN_FST", hue="Pop",order=['Model1',"Model10","Model13","Model4","Model7"],hue_order=['2pop','3pop'])
violin
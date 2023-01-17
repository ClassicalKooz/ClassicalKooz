import seaborn
import pandas as pd
df=pd.read_csv('/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/All_Models_PI_combined')
violin = seaborn.violinplot(data=df, x="Model",y="PI", hue="Pop",order=['Model1',"Model2","Model3","Model4","Model5","Model6"],hue_order=['2pop','3pop'])
violin
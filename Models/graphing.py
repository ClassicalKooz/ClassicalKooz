import seaborn
import pandas as pd
df=pd.read_csv('/home/michael/GhostPopulationSimulation/Models/All_Models_Fst_combined')
violin = seaborn.violinplot(data=df, x="Model",y="MEAN_FST", hue="Pop",order=['Model1',"Model2","Model3","Model4","Model5","Model6","Model7","Model8","Model9","Model10","Model11","Model12","Model13","Model14","Model15"],hue_order=['2pop','3pop'])
violin
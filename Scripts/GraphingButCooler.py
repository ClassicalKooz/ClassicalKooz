import seaborn 
import pandas as pd
from matplotlib import pyplot as plt
df=pd.read_csv('/home/michael/GhostPopulationSimulation/Models/All_Models_Pi_combined')
#violin = seaborn.violinplot(data=df, x="Model",y="PI", hue="Pop",order=['Model1',"Model2","Model3","Model4","Model5","Model6","Model7","Model8","Model9","Model10","Model11","Model12","Model13","Model14","Model15"],hue_order=['2pop','3pop'])
fig, axes = plt.subplots(3, 3, figsize=(18, 10))

fig.suptitle('Bilateral Models')
seaborn.violinplot(ax=axes[0, 0], data=df.loc[df["Model"] == 'Model1' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[0, 1], data=df.loc[df["Model"] == 'Model2' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[0, 2], data=df.loc[df["Model"] == 'Model3' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 0], data=df.loc[df["Model"] == 'Model4' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 1], data=df.loc[df["Model"] == 'Model5' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 2], data=df.loc[df["Model"] == 'Model6' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[2, 0], data=df.loc[df["Model"] == 'Model9' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[2, 1], data=df.loc[df["Model"] == 'Model8' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[2, 2], data=df.loc[df["Model"] == 'Model9' ], x='Pop', y='PI', order=['2pop','3pop'])

fig, axes = plt.subplots(2, 3, figsize=(18, 10))

fig.suptitle('Unilateral Models')
seaborn.violinplot(ax=axes[0, 0], data=df.loc[df["Model"] == 'Model10' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[0, 1], data=df.loc[df["Model"] == 'Model11' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[0, 2], data=df.loc[df["Model"] == 'Model12' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 0], data=df.loc[df["Model"] == 'Model13' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 1], data=df.loc[df["Model"] == 'Model14' ], x='Pop', y='PI', order=['2pop','3pop'])
seaborn.violinplot(ax=axes[1, 2], data=df.loc[df["Model"] == 'Model15' ], x='Pop', y='PI', order=['2pop','3pop'])
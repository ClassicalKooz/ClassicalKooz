import os
import pandas as pd
import glob

        
fst_files = glob.glob('Melissa_Simulated_Data/*PI_*_combined')
print(fst_files)
df_concat = pd.concat([pd.read_csv(f) for f in fst_files], ignore_index=True)
print(df_concat)
df_concat.to_csv("/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/All_Models_PI_combined",index=False)

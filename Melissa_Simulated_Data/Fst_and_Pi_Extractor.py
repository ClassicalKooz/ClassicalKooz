import os
import pandas as pd
import glob
for i in range (1,7):
    for j in range (2,4):
        
        fst_files = glob.glob('Melissa_Simulated_Data/Melissa_Model%d/PI/%dpop/*.pi'%(i,j))
        #print(fst_files)
        df_concat = pd.concat([pd.read_csv(f, sep="\t") for f in fst_files], ignore_index=True)
        print(df_concat)
        df_concat.to_csv("/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/Model%d_PI_%dpop_combined"%(i,j), index=False)

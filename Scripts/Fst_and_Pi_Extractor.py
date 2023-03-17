import os
import pandas as pd
import glob
for j in range (1,16):
    for i in range (2,4):      
        fst_files = glob.glob('Models/Model%d/PI/*.%dpop.*'%(j,i))
        #print(fst_files)
        df_concat = pd.concat([pd.read_csv(f, sep='\t') for f in fst_files], ignore_index=True)
        #print(df_concat)
        df_concat.insert(0,'Model','Model%d'%(j))
        df_concat.insert(1,'Pop','%dpop'%(i))
        print(df_concat)
        df_concat.to_csv("/home/michael/GhostPopulationSimulation/Models/Model%d_PI_%dpop_combined"%(j,i),index=False)

pi_files=glob.glob('Models/Model*_PI_*')
df_concat = pd.concat([pd.read_csv(f) for f in pi_files], ignore_index=True)
print(df_concat)
df_concat.to_csv("/home/michael/GhostPopulationSimulation/Models/All_Models_Pi_combined",index=False)
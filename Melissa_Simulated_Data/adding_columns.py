import pandas as pd
for i in range (2,7):
    for j in range(2,4):
        df = pd.read_csv("/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/Model%d_Fst_%dpop_combined"%(i,j))
        df.insert(0,'Model','Model%d'%(i))
        df.insert(1,'Pop','%dpop'%(j))
        print(df)
        df.to_csv("/home/michael/GhostPopulationSimulation/Melissa_Simulated_Data/Model%d_Fst_%dpop_combined"%(i,j),index=False)
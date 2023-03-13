import pandas as pd
import os   
path="/home/michael/GhostPopulationSimulation/models2/model1/"
dir_list = os.listdir(path)
#print(dir_list)
for FILE in dir_list:
    if 'windowed.weir.fst' in FILE:
        if '.log' not in FILE:
            db=pd.read_csv(path+FILE, sep='\t')
            #data_top = db.head()
            #print(data_top)
            #for col in db.columns:
                #print(col)
            print(FILE)
            print (db['WEIGHTED_FST'])
            print (db['MEAN_FST']) 
            dd=pd.DataFrame(db['WEIGHTED_FST'],db['MEAN_FST'])
                
                    
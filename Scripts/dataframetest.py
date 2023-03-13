import numpy as np
import pandas as pd

A=[]
B=[]

A.append(3*2)
B.append(4*5)
d={"A":A,"B":B}
df=pd.DataFrame(data=d)
print(df)
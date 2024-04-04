import pandas as pd
import numpy as np


phase =[]

t0 = 51197.0
F0 = 19.80244383176
F1 = -1.878039597E-10
F2 = 3.752027E-21  

df = pd.read_csv("/Users/paulmoran/Desktop/Coding/PulsarCodes/PSRB0540/Dec30_B0540_ts_bivot_fobs.txt", sep='\s+', header=None)
      
t = df.iloc[:,0]

for i in range(len(df)):
	tn = (t0 - t[i])*86400.0

	phase1 = F0*tn + 0.5*tn*tn*F1 + (1.0/6.0)*tn*tn*tn*F2

	phase.append((phase1 - (floor(phase1) + 1))*-1.0)


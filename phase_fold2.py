import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

divisor = 10

df = pd.read_csv("/Users/paulmoran/Desktop/Coding/PulsarCodes/PSRB0540/Dec30_B0540_ts_bivot_fobs.txt", sep='\s+', header=None)

phase = df.iloc[:,0]	# Take the phase part of the time-series
counts = df.iloc[:,1]	# Take the flux part of the time-series

#counts = np.genfromtxt("/Users/paulmoran/Desktop/Coding/PulsarCodes/Crab/Power_Spectra/crab_isaac_12ms.txt", delimiter='\s+')
#phase = np.genfromtxt("/Users/paulmoran/Desktop/Coding/PulsarCodes/Crab/Timing/12ms_03-02-2012/03-02-2012_phase_12ms_bary_Jodrell-Bank_2.txt", delimiter='\s+')

lightcurve = np.zeros(divisor)

for i in range(len(phase)):
	if 0.0 <= phase[i] <= 0.1:
		lightcurve[0] += counts[i]
	elif 0.1 <= phase[i] <= 0.2:
		lightcurve[1] += counts[i]
	elif 0.2 <= phase[i] <= 0.3:
		lightcurve[2] += counts[i]
	elif 0.3 <= phase[i] <= 0.4:
		lightcurve[3] += counts[i]
	elif 0.4 <= phase[i] <= 0.5:
		lightcurve[4] += counts[i]
	elif 0.5 <= phase[i] <= 0.6:
		lightcurve[5] += counts[i]
	elif 0.6 <= phase[i] <= 0.7:
		lightcurve[6] += counts[i]
	elif 0.7 <= phase[i] <= 0.8:
		lightcurve[7] += counts[i]	
	elif 0.8 <= phase[i] <= 0.9:
		lightcurve[8] += counts[i]
	elif 0.9 <= phase[i] <= 1.0:
		lightcurve[9] += counts[i]

bins = np.arange(-0.5,2.5,1/divisor)

#print(bins)
lightcurve2 = np.hstack((lightcurve, lightcurve,lightcurve))

plt.step(bins,lightcurve2)
plt.xlabel('Phase')
plt.ylabel('Counts')
plt.xlim(0,2)
plt.show()


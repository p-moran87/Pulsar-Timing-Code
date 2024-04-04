import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Define number of bins to use
divisor = 10

df = pd.read_csv("/Users/paulmoran/Desktop/Coding/PulsarCodes/PSRB0540/Dec30_B0540_ts_bivot_fobs.txt", sep='\s+', header=None)
phase = df.iloc[:,0]	# Take the phase part of the time-series
counts = df.iloc[:,1]	# Take the flux part of the time-series

#counts = np.genfromtxt("/Users/paulmoran/Desktop/Coding/PulsarCodes/Crab/Power_Spectra/crab_isaac_12ms.txt", delimiter='\s+')
#phase = np.genfromtxt("/Users/paulmoran/Desktop/Coding/PulsarCodes/Crab/Timing/12ms_03-02-2012/03-02-2012_phase_12ms_bary_Jodrell-Bank_2.txt", delimiter='\s+')

lightcurve = stats.binned_statistic(phase,counts, 'sum', bins=divisor)

bins = np.arange(-0.45,2.55,1/divisor)

print(len(phase))
print(len(counts))

#print(lightcurve)
lightcurve2 = np.hstack((lightcurve[0], lightcurve[0],lightcurve[0]))

plt.step(bins,lightcurve2)
plt.xlabel('Phase')
plt.ylabel('Counts')
plt.xlim(0,2)
plt.savefig('PSRB0540_salt_bivot_30Dec_lc.png')
plt.show()


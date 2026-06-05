import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats

# Define number of bins to use
divisor = 10

# Read in ephemeris file content
df = pd.read_csv(r"C:\Users\PaulMoran\OneDrive - Atlantic TU\Desktop\Research\Pulsar_Codes\Crab_VLT_ISAAC\Crab_JodrellBank_E_55889.txt", sep='\s+', header=None)
t0 = df.iloc[:,1][0]
F0 = df.iloc[:,1][1]
F1 = df.iloc[:,1][2]
F2 = df.iloc[:,1][3]

# PSRB0540-69 SALT/BVIT
#t0 = 56657.429 # ephemeris
#F0 = 19.80244383176 #old freq
#F0 = 19.71005842300069 #(measured PS) 19.73823 Hz estimate
#F1 = -1.878039597E-10
#F2 = 3.752027E-21

df = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\Crab_VLT_ISAAC\\03-02-2012_crab_12ms_bary_data.txt", sep='\s+', header=None)
t = df.iloc[:,0]
counts = df.iloc[:,2]	# Take the flux part of the time-series

#df = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\PSRB0540_SALT\\Dec30_B0540_full-ts_fobs_bary.txt", sep='\s+', header=None)
#t = df.iloc[:3500000,0]
#counts  = df.iloc[:3500000,2]

t_sec = (t0 - t) * 86400.0
#t_sec -= np.mean(t_sec)

phi = (
    F0 * t_sec
    + 0.5 * F1 * t_sec**2
    + (1.0 / 6.0) * F2 * t_sec**3
)

#phase = np.mod(phi * 0.999733, 1.0)
phase = np.mod(phi, 1.0)

#with open('Crab_VLT_ISAAC_03Feb12_PSFreq_Phases.txt', 'w') as file:
#    for phase in phase:
#        file.write(f"{phase}\n")

lightcurve = stats.binned_statistic(phase,counts, 'sum', bins=divisor)

bins = np.arange(-0.25,2.75,1/divisor)

y_error = np.sqrt(sum(counts)/divisor)

lightcurve2 = np.hstack((lightcurve[0], lightcurve[0],lightcurve[0]))

plt.step(bins,lightcurve2-min(lightcurve2))
#plt.errorbar(0.3, lightcurve[0][2], yerr=y_error, fmt="-", color="b")
plt.xlabel('Phase')
plt.ylabel('Counts')
plt.xlim(0,2)
#plt.savefig('Crab_vlt_isaac_03Feb12_lc_new_10bins.png')
#plt.savefig('PSRB0540_salt_bivot_30Dec_lc_new_30bins.png')
plt.show()
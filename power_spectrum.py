### SCRIPT FOR TAKING THE POWER SPECTRUM OF A SIGNAL
### check for periodicity at correct frequency

import matplotlib.pyplot as plt
from scipy.fftpack import rfft, rfftfreq
import pandas as pd
import numpy as np

flux = np.genfromtxt("/Users/paulmoran/Desktop/Coding/PulsarCodes/Crab/Power_Spectra/crab_isaac_12ms.txt", delimiter='\s+')
#df = pd.read_csv("/Users/paulmoran/Desktop/Coding/PulsarCodes/PSRB0540/Dec30_B0540_ts_bivot_fobs.txt", sep='\s+', header=None)
#flux = df.iloc[:,1]	# Take the flux part of the time-series

sample_rate = 0.012

# Get the frequency part of signal (real part of FFT)
freq = rfftfreq(flux.size, d=sample_rate)

# Power spectrum is the square of the real part of FFT
p = rfft(flux) **2

plt.plot(freq, p)
plt.xlim(0,50)
plt.ylim(0,400000000000)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
plt.savefig('Crab_vlt_isaac_12ms_ps_03Feb12.png')
plt.show()

# Get the frequency associated with the max value of the power spectrum 
max_index = np.argmax(p[100:])

print("The frequency (Hz) is: {}".format(freq[max_index+100]))

### SCRIPT FOR TAKING THE POWER SPECTRUM OF A SIGNAL
### check for periodicity at correct frequency

import matplotlib.pyplot as plt
from scipy.fftpack import rfft, rfftfreq
import pandas as pd
import numpy as np
import statistics as s

df = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\PSRB0540_69_counts_30Dec.txt", sep="\s+", header=None)
flux = df.iloc[:,1]
#flux = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\Pulsar_Codes\\Crab_isaac_12ms_03Feb12.txt", sep="\s+", header=None)
#flux = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\Crab_VLT_ISAAC\Region_1_crab_fluxes.txt", sep="\s+", header=None)

#sample_rate = 1/1840 # ESO (GASP)
#sample_rate = 0.012 # VLT/ISAAC (12 ms)
sample_rate = 0.0005 # SALT BVIT

# Get the frequency part of signal (real part of FFT)
freq = rfftfreq(flux.size, d=sample_rate)

# Power spectrum is the square of the real part of FFT
power = rfft(flux) **2

plt.plot(freq, power)
plt.xlim(0,50)
plt.ylim(0,2e8)
#plt.ylim(0,4e11) # Crab VLT/ISAAC
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
#plt.title("Crab VLT/ISAAC f=29.71006 Hz, 700 sigma")
#plt.savefig('Crab_VLT_ISAAC_03Feb12_ps.png')
plt.show()

# Remove spurious signals
for i in range(len(power)):
    if power[i] >= 1e8:
        if 0.0 <= freq[i] <= 3.0:
            power[i] = s.mean(power[200:500])
        elif 10.0 <= freq[i] <= 18.0:
            power[i] = s.mean(power[200:500])
        elif 22.0 <= freq[i] <= 25.0:
            power[i] = s.mean(power[200:500])
        elif 31.0 <= freq[i] <= 35.0:
            power[i] = s.mean(power[200:500])
            
plt.plot(freq, power)
plt.xlim(0,50)
plt.ylim(0,2e8)  # PSRB0540 SALT
#plt.ylim(0,4e11) # Crab VLT/ISAAC
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
#plt.title("Crab VLT/ISAAC f=29.70312 Hz, 700 sigma")
#plt.savefig('Crab_VLT_ISAAC_03Feb12_ps_clean.png')
plt.show()


# Significance of signal: signal - mean / std dev
std_bkg = s.stdev(power[200:500])
mean_bkg = s.mean(power[200:500])
sigma =  (np.max(power) - mean_bkg) / std_bkg
sigma = round(sigma,2)

print("Frequency range:", freq.min(), "to", round(freq.max(),5), "Hz")
print("Frequency resolution:", round(freq.max()/len(freq),5), "Hz")
print("Peak frequency:", freq[np.argmax(power)], "Hz")
print("Significance", sigma)







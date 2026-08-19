# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 16:04:37 2026

@author: PaulMoran
"""

### Determine the Fourier Power Spectrum of a transient source
### Filter out alias frequencies

from pathlib import Path

import matplotlib.pyplot as plt
from scipy.fftpack import rfft, rfftfreq
import pandas as pd
import numpy as np
import statistics as s

path = input("Enter file path: (e.g. 'C:/Users/time-series.txt') \n").strip().strip('"')

df = pd.read_csv(path, sep='\s+', header=None)

flux = df.iloc[:,2]
flux = flux.values

signal = flux - np.mean(flux)

sample_rate = float(input("Please enter frame rate of the instrument in seconds:\n"))
mains_freq = float(input("Please enter frequency of the mains current in Hz (e.g. 50 or 60 Hz:)\n"))

target_freq = float(input("Please enter frequency of transient target (rotation or period) Hz:)\n"))
target_name = input("Please enter the name of the astronomical target?:\n")

sample_freq = 1.0 / sample_rate # sampling frequency

# Get the frequency part of signal (real part of FFT)
freq = rfftfreq(signal.size, d=sample_rate)

# Power spectrum is the square of the real part of FFT
power = np.abs(rfft(signal)) **2

# Determine limits for x and y axis of FFT PS
x_max = freq.max()
y_max = power.max()*1.1

plt.plot(freq, power)
plt.xlim(0,x_max)
plt.ylim(0,y_max) 
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
plt.savefig('Raw_PS_{target_name}.png')
plt.show()

# Remove spurious signals
for i in range(len(power)):
    if power[i] >= 1e8:
        if 0.0 <= freq[i] <= 3.0:
            power[i] = s.mean(power[200:2000])
    for k in range(1,10):
        f_alias_mains = np.abs(sample_freq - k*mains_freq) 
        if (f_alias_mains - 3.0) <= freq[i] <= (f_alias_mains + 3.0):
            power[i] = s.mean(power[200:2000])
        f_alias_target = np.abs(sample_freq - k * target_freq)
        if (f_alias_target - 3.0) <= freq[i] <= (f_alias_target + 3.0):
            power[i] = s.mean(power[200:2000])
            
plt.plot(freq, power)
plt.xlim(0,x_max)
plt.ylim(0,y_max) 
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
plt.savefig('Clean_PS_{target_name}.png')
plt.show()

# Significance of signal: signal - mean / std dev
std_bkg = s.stdev(power[200:2000])
mean_bkg = s.mean(power[200:2000])
sigma =  (np.max(power) - mean_bkg) / std_bkg
sigma = round(sigma,2)

print("Frequency range:", freq.min(), "to", round(freq.max(),5), "Hz")
print("Frequency resolution:", round(2*freq.max()/len(freq),5), "Hz")
print("Peak frequency:", freq[np.argmax(power)], "Hz")
print("Significance", sigma)







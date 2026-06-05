# -*- coding: utf-8 -*-
"""
Created on Wed May 27 12:53:42 2026

@author: PaulMoran
"""
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation
from numpy.fft import rfft, rfftfreq
import matplotlib.pyplot as plt
import pandas as pd
import statistics as s

# USER INPUTS
# ---------------------------
df = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\Crab_VLT_ISAAC\\03-02-2012_crab_12ms_bary_data.txt", sep='\s+', header=None)
flux = df.iloc[:,2]	# Take the flux part of the time-series

#df = pd.read_csv("C:\\Users\\PaulMoran\\OneDrive - Atlantic TU\\Desktop\\Research\\Pulsar_Codes\\PSRB0540_69_counts_30Dec.txt", sep="\s+", header=None)
#flux = df.iloc[:,1]

#sample_rate = 1/1840 # ESO (GASP)
sample_rate = 0.012 # VLT/ISAAC (12 ms)
#sample_rate = 0.0005 # SALT BVIT

t0_mjd = 55960.093826  # start time of VLT/ISAAC observations (MJD)
#t0_mjd = 56657.429  # start time of SALT/BVIT observations (MJD)

# Observatory (example: Paranal / VLT — change if needed)
vlt = EarthLocation.of_site('paranal')

salt = EarthLocation(
    lat=-32.3794 * u.deg,
    lon=20.8107 * u.deg,
    height=1798 * u.m
)

# Crab pulsar coordinates
crab = SkyCoord(ra=83.6331 * u.deg, 
                    dec=22.0145 * u.deg,frame='icrs')

psr_b0540 = SkyCoord(ra=85.0467 * u.deg,
                     dec=-69.3319 * u.deg, frame="icrs")

# STEP 1: Create time array
N = flux.size
times_sec = (np.arange(N) + 0.5) * sample_rate
times_mjd = t0_mjd + times_sec / 86400.0

# STEP 2: Apply barycentric correction
t = Time(times_mjd, format='mjd', scale='utc', location=vlt)

# Light travel time correction
ltt_bary = t.light_travel_time(crab)

# Convert to barycentric dynamical time (TDB)
t_bary = t.tdb + ltt_bary

# STEP 3: Convert to seconds relative to start
dt_bary = (t_bary - t_bary[0]).to(u.s).value

# STEP 4: Compute FFT in barycentric frame
# Use mean spacing (nearly uniform after correction)
dt_mean = np.mean(np.diff(dt_bary))

# Get the frequency part of signal (real part of FFT)
freq = rfftfreq(flux.size, d=sample_rate)

# Power spectrum is the square of the real part of FFT
power = np.abs(rfft(flux)) **2

nyquist = int(freq.max())
ymax = np.max(power[100:])

# OUTPUT
plt.plot(freq, power)
plt.xlim(0,nyquist)
plt.ylim(0,5e11)
plt.xlabel("Frequency (Hz)")
plt.ylabel("Power")
#plt.title("Crab VLT/ISAAC f=29.71006 Hz, 700 sigma")
#plt.savefig('Crab_VLT_ISAAC_03Feb12_ps.png')
plt.show()

#peaks = []
#for i in range(len(power)):
#    if power[i] > 1e10:
#        peaks.append(freq[i])
#print(peaks)

# Remove spurious signals
for i in range(len(power)):
    if power[i] >= 1e10:
        if 0.0 <= freq[i] <= 3.0:
            power[i] = s.mean(power[200:500])
        elif 10.0 <= freq[i] <= 18.0:
            power[i] = s.mean(power[200:500])
        elif 22.0 <= freq[i] <= 25.0:
            power[i] = s.mean(power[200:500])
        elif 31.0 <= freq[i] <= 38.0:
            power[i] = s.mean(power[200:500])
            
plt.plot(freq, power)
plt.xlim(0,nyquist)
#plt.ylim(0,2.5e8)
plt.ylim(0,ymax) # Crab VLT/ISAAC
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

# Read in ephemeris file content
df = pd.read_csv(r"C:\Users\PaulMoran\OneDrive - Atlantic TU\Desktop\Research\Pulsar_Codes\Crab_VLT_ISAAC\Crab_JodrellBank_E_55889.txt", sep='\s+', header=None)
t_ephem = df.iloc[:,1][0]
F0 = df.iloc[:,1][1]
F1 = df.iloc[:,1][2]
F2 = df.iloc[:,1][3]

dt_sec = (t_bary.mjd - t_ephem) * 86400.0
#dt_sec = dt_sec - np.mean(dt_sec)

# Frequency search around expected Crab value
trial_freqs = np.linspace(29.702, 29.705, 5000)  #Crab
#trial_freqs = np.linspace(19.73, 19.75, 5000)  #PSRB0540

chi2_vals = []

for f in trial_freqs:
    phi = (
        f * dt_sec
        + 0.5 * F1 * dt_sec**2
        + (1/6) * F2 * dt_sec**3
    )

    phase = np.mod(phi, 1.0)

    hist, _ = np.histogram(phase, bins=50)

    mean = np.mean(hist)
    chi2 = np.sum((hist - mean)**2 / mean)
    chi2_vals.append(chi2)

best_freq = trial_freqs[np.argmax(chi2_vals)]

# sub-grid precision
i = np.argmax(chi2_vals)

x = trial_freqs[i-1:i+2]
y = np.array(chi2_vals[i-1:i+2])

coeffs = np.polyfit(x, y, 2)
f_peak = -coeffs[1] / (2 * coeffs[0])

print("Frequency range:", freq.min(), "to", round(freq.max(),5), "Hz")
print("Frequency resolution:", round(freq.max()/len(freq),5), "Hz")
print("Peak frequency:", freq[np.argmax(power)], "Hz")
print("Significance", sigma)
print("Best frequency:", best_freq, "Hz")
print("Refined frequency:", f_peak, "Hz")



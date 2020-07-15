
"""

Data: Temperature and Salinity time series from SIO Scripps Pier
    Salinity: measured in PSU at the surface (~0.5m) and at depth (~5m)
    Temp: measured in degrees C at the surface (~0.5m) and at depth (~5m)
- Timestamp included beginning in 1990

"""

# imports
import sys,os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime 
from scipy import signal
import scipy.stats as ss
import SIO_modules as SIO_mod
from importlib import reload
reload(SIO_mod)

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
ENSO_data = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_data.xlsx')
ENSO_data_recent = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_recent_data.xlsx')
PDO_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_PDO_data.csv', skiprows = 1)

# path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

# convert year, month, day columns to single DATE column
sal_data['DATE'] = pd.to_datetime(sal_data[['YEAR', 'MONTH', 'DAY']])
temp_data['DATE'] = pd.to_datetime(temp_data[['YEAR', 'MONTH', 'DAY']])
ENSO_data_all = ENSO_data.append(ENSO_data_recent[323:], ignore_index = True)
PDO_data['DATE'] = pd.to_datetime(PDO_data['Date'], format='%Y%m')

# remove uncertain data(SURF_FLAG between 1 and 4), replace with NaN, then interpolate 
for i in range(0,len(sal_data['SURF_SAL_PSU'])):
	if (sal_data['SURF_FLAG'][i] >= 1) and (sal_data['SURF_FLAG'][i] <=4):
		sal_data['SURF_SAL_PSU'][i] = np.nan

for i in range(0,len(temp_data['SURF_TEMP_C'])):
	if (sal_data['SURF_FLAG'][i] >= 1) and (sal_data['SURF_FLAG'][i] <=4):
		sal_data['SURF_SAL_PSU'][i] = np.nan

# interpolate missing temp and sal data
sal_data['SURF_SAL_PSU'] = sal_data['SURF_SAL_PSU'].interpolate()
temp_data['SURF_TEMP_C'] = temp_data['SURF_TEMP_C'].interpolate()
sal_data['SURF_SAL_PSU'][0] = sal_data['SURF_SAL_PSU'][1]

# remove the average from the sal and temp data and create new columns
sal_data['SURF_SAL_PSU_NOAVG'] = sal_data['SURF_SAL_PSU'] - sal_data['SURF_SAL_PSU'].mean()
temp_data['SURF_TEMP_C_NOAVG'] = temp_data['SURF_TEMP_C'] - temp_data['SURF_TEMP_C'].mean()

# remove trends from the sal and temp data and create new columns
sal_fit = np.polyfit(sal_data.index,sal_data['SURF_SAL_PSU_NOAVG'],1)
sal_fit_fn = np.poly1d(sal_fit)
temp_fit = np.polyfit(temp_data.index,temp_data['SURF_TEMP_C_NOAVG'],1)
temp_fit_fn = np.poly1d(temp_fit)

sal_fit_value = sal_fit_fn(sal_data.index)
temp_fit_value = temp_fit_fn(temp_data.index)

sal_data['SURF_SAL_PSU_DETREND'] = sal_data['SURF_SAL_PSU_NOAVG'] - sal_fit_value
temp_data['SURF_TEMP_C_DETREND'] = temp_data['SURF_TEMP_C_NOAVG'] - temp_fit_value

# butterworth low pass filter for temperature and salinity
fs = 1 # sampling frequency, once per day
fc = 1/500 # cut-off frequency of the filter (cut off periods shorter than 500 days)
w = fc / (fs / 2) #normalize the frequency
b, a = signal.butter(4, w, 'low')
temp_output = signal.filtfilt(b, a, temp_data['SURF_TEMP_C_DETREND'])
sal_output = signal.filtfilt(b, a, sal_data['SURF_SAL_PSU_DETREND'])

# create dataframe with spectra for each variable
spectra_temp_df = pd.DataFrame(columns = ['Temp_freq', 'Temp_spec', 'Temp_fft'])
spectra_sal_df = pd.DataFrame(columns = ['Sal_freq', 'Sal_spec', 'Sal_fft'])
spectra_PDO_df = pd.DataFrame(columns = ['PDO_freq', 'PDO_spec', 'PDO_fft'])
spectra_ENSO_df = pd.DataFrame(columns = ['ENSO_freq', 'ENSO_spec', 'ENSO_fft'])

# for coherence, start all records at 1916-01-01
# ENSO data [20:] 1916-09-01 onward, monthly// ends now, through 2019-05-01 [:1254]
# Temp data [10:] 1916-09-01 onward, daily // ends 2019-05-31
# PDO data [752:] 1916-09-01 onward, monthly// ends now, thorugh 2019-05-01 [:1985]

# compute spectral variables for each variable
for j in range(0,4):
	data_sets = [temp_data['SURF_TEMP_C_DETREND'][10:], sal_data['SURF_SAL_PSU_DETREND'], PDO_data['Value'][744:], ENSO_data_all['VALUE'][20:1255]]
	freq, spec, spec_amp, fft, delt, freq_T, freq_nyquist = SIO_mod.var_fft(data_sets[j])
	if j == 0:
		spectra_temp_df['Temp_freq'] = freq
		spectra_temp_df['Temp_spec'] = spec
		spectra_temp_df['Temp_fft'] = fft
	if j == 1:
		spectra_sal_df['Sal_freq'] = freq
		spectra_sal_df['Sal_spec'] = spec
		spectra_sal_df['Sal_fft'] = fft
	if j == 2:
		spectra_PDO_df['PDO_freq'] = freq
		spectra_PDO_df['PDO_spec'] = spec
		spectra_PDO_df['PDO_fft'] = fft
	if j == 3:
		spectra_ENSO_df['ENSO_freq'] = freq
		spectra_ENSO_df['ENSO_spec'] = spec
		spectra_ENSO_df['ENSO_fft'] = fft
delt = 1
n_av = 30
# temp_spec_amp_av, temp_spec_phase_av, temp_freq_av, temp_count = SIO_mod.band_average(spectra_temp_df['Temp_fft'][0:625], spectra_temp_df['Temp_fft'][0:625], spectra_temp_df['Temp_freq'][0:6], n_av, delt)
n_av = 1
# cospec_amp,cospec_phase,freq_av,count=SIO_mod.band_average(spectra_temp_df['Temp_fft'], spectra_ENSO_df['ENSO_fft'], spectra_ENSO_df['ENSO_freq'], n_av, delt)
# coh_sq=cospec_amp**2/(temp_spec*sal_spec)
tstr = 'Temp Spectra and Temp Spectra Band Av = 30'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].loglog(temp_freq_av, temp_spec_amp_av)
axes[0].set_xlabel('$\omega$ (radians/day)')
axes[0].set_ylabel('Temp Band averaged Spec n_av = 30')
axes[0].set_ylim(10**-8, 10**5)
axes[1].set_xlabel('$\omega$ (radians/day)')
axes[1].set_ylabel('Temp Spec no band average')
axes[1].loglog(spectra_temp_df['Temp_freq'], spectra_temp_df['Temp_spec'])
axes[1].loglog(spectra_ENSO_df['ENSO_freq'], spectra_ENSO_df['ENSO_spec'], color = 'red')
axes[1].set_ylim(10**-8, 10**5)
fig.suptitle(tstr)
plt.show()

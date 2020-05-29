"""

Data: Temeprature and Salinity time series from SIO Scripps Pier
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

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
ENSO_data = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_data.xlsx')
ENSO_data_recent = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_recent_data.xlsx')
precip_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_Precip_data.csv')
PDO_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_PDO_data.csv', skiprows = 1)

# path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

# convert year, month, day columns to single DATE column
sal_data['DATE'] = pd.to_datetime(sal_data[['YEAR', 'MONTH', 'DAY']])
temp_data['DATE'] = pd.to_datetime(temp_data[['YEAR', 'MONTH', 'DAY']])
precip_data['Date'] = pd.to_datetime(precip_data['DATE'], format='%Y-%m-%d')
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

# filter data
sal_ma = sal_data['SURF_SAL_PSU_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'boxcar').mean()
sal_tri = sal_data['SURF_SAL_PSU_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'triang').mean()
temp_ma = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'boxcar').mean()
temp_tri = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'triang').mean()

# butterworth low pass filter for temperature and salinity
fs = 1 # sampling frequency, once per day
fc = 1/500 # cut-off frequency of the filter (cut off periods shorter than 500 days)
w = fc / (fs / 2) #normalize the frequency
b, a = signal.butter(4, w, 'low')
temp_output = signal.filtfilt(b, a, temp_data['SURF_TEMP_C_DETREND'])
sal_output = signal.filtfilt(b, a, sal_data['SURF_SAL_PSU_DETREND'])

# create variables to carry out spectral analysis (SIO temp, SIO sal, PDO, ENSO)
def var_fft(j):
	data_sets = [temp_data['SURF_TEMP_C_DETREND'], sal_data['SURF_SAL_PSU_DETREND'], PDO_data['Value'][744:], ENSO_data_all['VALUE']]
	ll = len(data_sets[j])
	if j <= 1:
		delt = 1
	if j > 1:
		delt = 30
	data_fft = data_sets[j]

	ll_half = int(ll/2+1)
	freq = np.zeros(ll_half)
	T_length = ll
	freq_nyquist = np.pi/delt
	freq_T = 2.*np.pi/(ll*delt)

	omega0 = 2.*np.pi/(T_length*delt)
	for i in range(0,ll_half):
		freq[i] = i*omega0

	fft = np.fft.rfft(data_fft,n=ll)
	fourier_amp = np.sqrt((np.real(fft)**2+np.imag(fft)**2))
	fourier_phase = 180.*np.arctan2(np.imag(fft),np.real(fft))/np.pi
	spec = np.real(fft*np.conj(fft))/(2.*np.pi*T_length*delt)
	spec_amp  = (np.absolute(fft))**2/(2*np.pi*T_length*delt)

	return(freq,spec,spec_amp,delt,freq_T,freq_nyquist)

# dictionaries to define legend names, colors, and labels
legend_name = ['Temp', 'Sal', 'PDO', 'ENSO']
color_sets = ['mediumaquamarine','cornflowerblue','green','blue']
spectra_label = ['($^\circ$C$^2$)/(cycles/day)', '(PSU$^2$)/(cycles/day)','(PDO Index$^2$)/(cycles/month)','(ENSO Index$^2$)/(cycles/month)']
title_label = ['Power Spectrum SIO Temeprature', 'Power Spectrum SIO Salinity', 'Power Spectrum PDO Index','Power Spectrum ENSO Index']
x_lab = ['Day', 'Day', 'Month', 'Month']

# plot spectra
plt.close()
for j in range(0,4):
	fig = plt.figure(figsize = (10,6))
	ax = fig.add_subplot(111)
	freq, spec, spec_amp, delt, freq_T, freq_nyquist = var_fft(j)
	ax.loglog(freq, spec, label = legend_name[j], color = color_sets[j])
	if j==1:
		ax.set_ylim(10**-8, 10**2)
	if j==0:
		ax.set_ylim(10**-7, 10**5)
	if j>=2:
		ax.autoscale()
	ax.axvline(freq_nyquist, color = 'black', linestyle = '--', alpha = 0.5)
	ax.text(0.9, 0.9,'$\omega_{max}$', alpha = 0.5, transform = ax.transAxes)
	ax.axvline(freq_T, color = 'black', linestyle = '--', alpha = 0.5)
	ax.text(0.05, 0.9,'$\omega_o$', alpha = 0.5, transform = ax.transAxes)
	ax.set_title(title_label[j])
	ax.set_ylabel('Energy Density \n' + spectra_label[j])
	ax.set_xlabel('Cycles per ' + x_lab[j])
	plt.show()



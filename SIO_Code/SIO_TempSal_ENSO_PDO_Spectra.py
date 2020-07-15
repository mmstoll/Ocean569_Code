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

# butterworth low pass filter for temperature and salinity
fs = 1 # sampling frequency, once per day
fc = 1/500 # cut-off frequency of the filter (cut off periods shorter than 500 days)
w = fc / (fs / 2) #normalize the frequency
b, a = signal.butter(4, w, 'low')
temp_output = signal.filtfilt(b, a, temp_data['SURF_TEMP_C_DETREND'])
sal_output = signal.filtfilt(b, a, sal_data['SURF_SAL_PSU_DETREND'])

# create variables to carry out spectral analysis (SIO temp and sal)
ssl = len(sal_data['SURF_SAL_PSU_NOAVG'])
tt = len(temp_data['SURF_TEMP_C_NOAVG'])
pi = np.pi
ss_half = int(ssl/2+1)
tt_half = int(tt/2+1)
ss_freq = np.zeros(ss_half)
tt_freq = np.zeros(tt_half)

ss_T_length = ssl
ss_delt = 1 #days
tt_T_length = tt
tt_delt = 1 #days

ss_freq_nyquist=pi/ss_delt
ss_freq_T=2.*pi/(ssl*ss_delt)
tt_freq_nyquist=pi/tt_delt
tt_freq_T=2.*pi/(tt*tt_delt)  
freq_ann = 2*pi/365

ss_omega0 = 2.*pi/(ss_T_length*ss_delt)
tt_omega0 = 2.*pi/(tt_T_length*tt_delt)
for i in range(0,ss_half):
    ss_freq[i]=i*ss_omega0
for i in range(0,tt_half):
	tt_freq[i] = i*tt_omega0

# create variables to carry out spectral analysis (ENSO and PDO)
ee = len(ENSO_data_all['VALUE'])
pp = len(PDO_data['Value'][744:])
ee_half = int(ee/2+1)
pp_half = int(pp/2+1)
ee_freq = np.zeros(ee_half)
pp_freq = np.zeros(pp_half)

ee_T_length = ee
pp_T_length = pp
ee_delt = 30 #days (measured monthly)
pp_delt = 30 #days (measured monthly)

ee_freq_nyquist=pi/ee_delt
ee_freq_T=2.*pi/(ee*ee_delt)
pp_freq_nyquist=pi/pp_delt
pp_freq_T=2.*pi/(pp*pp_delt)  

ee_omega0 = 2.*pi/(ee_T_length*ee_delt)
pp_omega0 = 2.*pi/(pp_T_length*pp_delt)
for i in range(0,ee_half):
    ee_freq[i]=i*ee_omega0
for i in range(0,pp_half):
	pp_freq[i] = i*pp_omega0

# compute the power spectrum for temperature and salinity, detrended
ss_fft = np.fft.rfft(sal_output,n=ssl)
ss_fourier_amp = np.sqrt((np.real(ss_fft)**2+np.imag(ss_fft)**2))
ss_fourier_phase=180.*np.arctan2(np.imag(ss_fft),np.real(ss_fft))/pi
ss_spec = np.real(ss_fft*np.conj(ss_fft))/(2.*pi*ss_T_length*ss_delt)
ss_spec_amp = (np.absolute(ss_fft))**2/(2*pi*ss_T_length*ss_delt)

tt_fft = np.fft.rfft(temp_data['SURF_TEMP_C_DETREND'],n=tt)
tt_fourier_amp = np.sqrt((np.real(tt_fft)**2+np.imag(tt_fft)**2))
tt_fourier_phase=180.*np.arctan2(np.imag(tt_fft),np.real(tt_fft))/pi
tt_spec = np.real(tt_fft*np.conj(tt_fft))/(2.*pi*tt_T_length*tt_delt)
tt_spec_amp = (np.absolute(tt_fft))**2/(2*pi*tt_T_length*tt_delt)

tt_fft_sm = np.fft.rfft(temp_output,n=tt)
tt_fourier_amp_sm = np.sqrt((np.real(tt_fft_sm)**2+np.imag(tt_fft_sm)**2))
tt_fourier_phase_sm=180.*np.arctan2(np.imag(tt_fft_sm),np.real(tt_fft_sm))/pi
tt_spec_sm = np.real(tt_fft_sm*np.conj(tt_fft_sm))/(2.*pi*tt_T_length*tt_delt)
tt_spec_amp_sm = (np.absolute(tt_fft_sm))**2/(2*pi*tt_T_length*tt_delt)

ee_fft = np.fft.rfft(ENSO_data_all['VALUE'],n=ee)
ee_fourier_amp = np.sqrt((np.real(ee_fft)**2+np.imag(ee_fft)**2))
ee_fourier_phase=180.*np.arctan2(np.imag(ee_fft),np.real(ee_fft))/pi
ee_spec = np.real(ee_fft*np.conj(ee_fft))/(2.*pi*ee_T_length*ee_delt)
ee_spec_amp = (np.absolute(ee_fft))**2/(2*pi*ee_T_length*ee_delt)

pp_fft = np.fft.rfft(PDO_data['Value'][744:],n=pp)
pp_fourier_amp = np.sqrt((np.real(pp_fft)**2+np.imag(pp_fft)**2))
pp_fourier_phase=180.*np.arctan2(np.imag(pp_fft),np.real(pp_fft))/pi
pp_spec = np.real(pp_fft*np.conj(pp_fft))/(2.*pi*pp_T_length*pp_delt)
pp_spec_amp = (np.absolute(pp_fft))**2/(2*pi*pp_T_length*pp_delt)

# find confidence intervals
# df = 10
# conf_l = np.zeros(len(tt_spec))
# conf_h = np.zeros(len(tt_spec))
# for i in range(0,len(tt_spec)):
# 	conf_l[i] = tt_spec[i] * df / ss.chi2.ppf([0.975], df)
# 	conf_h[i] = tt_spec[i] * df / ss.chi2.ppf([0.025], df)

# conf = np.zeros(len(tt_spec))
# for i in range (0,len(tt_spec)):
# 	conf[i] = conf_h[i] - conf_l[i]
# df = 10
# conf_lim = 0.95
# conf_above = (1.-conf_lim)/2
# conf_below = 1.-(1.-conf_lim)/2
# mark_rt = stats.chi2.ppf(conf_below,df)
# mark_lf = stats.chi2.ppf(conf_above,df)
# yy_upper_nolog = mark_rt
# yy_lower_nolog = mark_lf
# yy_upper0[i] = (yy_upper_nolog - df)/df
# yy_lower0[i] = (df - yy_lower_nolog)/df

# fig = plt.figure()
# plt.loglog(tt_freq,conf, color = 'k')
# plt.loglog(tt_freq, conf_l, color = 'blue')
# plt.loglog(tt_freq, conf_h, color = 'red')
# plt.show()

# plot the spectra, detrended
tstr = 'SIO Temperature Power Spectra, ENSO and PDO'
im_name = 'SIO_Temp_ENSOPDO_Spectra.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].loglog(tt_freq,tt_spec_amp, color='mediumaquamarine')
axes[0].loglog(ee_freq,ee_spec_amp, color = 'black')
# axes[0].fill_between(tt_freq, conf_l, conf_h, color='k', alpha = 0.2)
axes[0].plot([tt_freq_nyquist,tt_freq_nyquist], [10**-8,10**6],'--k', alpha = 0.5)
axes[0].text(10**0.16,10**4,'$\omega_{max}$', alpha = 0.5)
axes[0].plot([tt_freq_T,tt_freq_T], [10**-8,10**6], '--k', alpha = 0.5)
axes[0].text(10**-3.7,10**4,'$\omega_o$', alpha = 0.5)
axes[0].plot([freq_ann,freq_ann], [10**-8,10**6],'--k', alpha = 0.5)
axes[0].set_ylim(10**-7, 10**5)
axes[0].set_xlabel('$\omega$ (radians/day)')
axes[0].set_ylabel('Energy Density \n($^\circ$C$^2$)/(radians/day)')
axes[1].loglog(tt_freq,tt_spec_amp_sm, color='mediumaquamarine')
axes[1].loglog(pp_freq,pp_spec_amp, color = 'black')
axes[1].plot([tt_freq_nyquist,tt_freq_nyquist], [10**-10,10**2], '--k', alpha = 0.5)
axes[1].text(10**0.16,10**4,'$\omega_{max}$', alpha = 0.5)
axes[1].plot([tt_freq_T,tt_freq_T], [10**-10,10**2], '--k', alpha = 0.5)
axes[1].text(10**-3.7,10**4,'$\omega_o$', alpha = 0.5)
axes[1].plot([freq_ann,freq_ann], [10**-10,10**2],'--k', alpha = 0.5)
axes[1].set_ylim(10**-7, 10**5)
axes[1].set_xlabel('$\omega$ (radians/day)')
axes[1].set_ylabel('Energy Density \n(PSU$^2$)/(radians/day)')
fig.suptitle(tstr)
# plt.savefig(path_out + im_name)
plt.show()




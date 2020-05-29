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

# create dataframe with spectra for each variable
spectra_temp_df = pd.DataFrame(columns = ['Temp_freq', 'Temp_spec'])
spectra_sal_df = pd.DataFrame(columns = ['Sal_freq', 'Sal_spec'])
spectra_PDO_df = pd.DataFrame(columns = ['PDO_freq', 'PDO_spec'])
spectra_ENSO_df = pd.DataFrame(columns = ['ENSO_freq', 'ENSO_spec'])

for j in range(0,4):
	freq, spec, spec_amp, delt, freq_T, freq_nyquist = var_fft(j)
	if j == 0:
		spectra_temp_df['Temp_freq'] = freq
		spectra_temp_df['Temp_spec'] = spec
	if j == 1:
		spectra_sal_df['Sal_freq'] = freq
		spectra_sal_df['Sal_spec'] = spec
	if j == 2:
		spectra_PDO_df['PDO_freq'] = freq
		spectra_PDO_df['PDO_spec'] = spec
	if j == 3:
		spectra_ENSO_df['ENSO_freq'] = freq
		spectra_ENSO_df['ENSO_spec'] = spec

# find confidence intervals
spec_set = [spectra_temp_df['Temp_spec'], spectra_sal_df['Sal_spec'], spectra_PDO_df['PDO_spec'], spectra_ENSO_df['ENSO_spec']]
freq_set = [spectra_temp_df['Temp_freq'], spectra_sal_df['Sal_freq'], spectra_PDO_df['PDO_freq'], spectra_ENSO_df['ENSO_freq']]

conf_temp_l = np.zeros(len(spec_set[0]))
conf_temp_h = np.zeros(len(spec_set[0]))
conf_sal_l = np.zeros(len(spec_set[1]))
conf_sal_h = np.zeros(len(spec_set[1]))
conf_PDO_l = np.zeros(len(spec_set[2]))
conf_PDO_h = np.zeros(len(spec_set[2]))
conf_ENSO_l = np.zeros(len(spec_set[3]))
conf_ENSO_h = np.zeros(len(spec_set[3]))

df = 2
for j in range (0,4):
	if j == 0:
		for i in range(0,len(spec_set[j])):
			conf_temp_l[i] = spec_set[j][i] * df / ss.chi2.ppf([0.975], df)
			conf_temp_h[i] = spec_set[j][i] * df / ss.chi2.ppf([0.025], df)
	if j == 1:
		for i in range(0,len(spec_set[j])):
			conf_sal_l[i] = spec_set[j][i] * df / ss.chi2.ppf([0.975], df)
			conf_sal_h[i] = spec_set[j][i] * df / ss.chi2.ppf([0.025], df)
	if j == 2:
		for i in range(0,len(spec_set[j])):	
			conf_PDO_l[i] = spec_set[j][i] * df / ss.chi2.ppf([0.975], df)
			conf_PDO_h[i] = spec_set[j][i] * df / ss.chi2.ppf([0.025], df)
	if j == 3:

		for i in range(0,len(spec_set[j])):
			conf_ENSO_l[i] = spec_set[j][i] * df / ss.chi2.ppf([0.975], df)
			conf_ENSO_h[i] = spec_set[j][i] * df / ss.chi2.ppf([0.025], df)

conf_l_set = [conf_temp_l, conf_sal_l, conf_PDO_l, conf_ENSO_l]
conf_h_set = [conf_temp_h, conf_sal_h, conf_PDO_h, conf_ENSO_h]

for j in range(0,4):
	fig = plt.figure(figsize = (10,6))
	ax = fig.add_subplot(111)
	ax.loglog(freq_set[j], spec_set[j], label = legend_name[j], color = color_sets[j])
	ax.fill_between(freq_set[j], conf_l_set[j], conf_h_set[j] , color = 'k', alpha = 0.2)
	ax.set_title(title_label[j])
	ax.set_ylabel('Energy Density \n' + spectra_label[j])
	ax.set_xlabel('Cycles per ' + x_lab[j])
	if j <= 1:
		ax.set_ylim(10**-8, 10**5)
	plt.show()

# specific for temperature
freq,spec,spec_amp,delt,freq_T,freq_nyquist = var_fft(0)

# conf = np.zeros(len(spec))
nn = len(spec)
xx=np.zeros(nn)
yy_lower0=np.zeros(nn)
yy_upper0=np.zeros(nn)
frac=0.1

for i in range(0,nn-2):
    ii=frac*(i+2)
    nu=ii
    d_mean=nu
    conf_lim=0.95
    xx[i]=ii
    conf_above=(1.-conf_lim)/2
    conf_below=1.-(1.-conf_lim)/2.
    mark_rt = ss.chi2.ppf(conf_below,nu)
    mark_lf = ss.chi2.ppf(conf_above,nu)
    yy_upper_nolog=mark_rt
    yy_lower_nolog=mark_lf
    yy_upper0[i]=(yy_upper_nolog-d_mean)/d_mean
    yy_lower0[i]=(d_mean-yy_lower_nolog)/d_mean
ind_set=nn-2
xx=xx[0:ind_set]
yy_lower1=yy_lower0[0:ind_set]
yy_upper1=yy_upper0[0:ind_set]

fig = plt.figure()
plt.plot(xx,yy_upper1)
plt.plot(xx, -yy_lower1)
plt.show()

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

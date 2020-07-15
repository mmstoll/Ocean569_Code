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
import SIO_modules as SIO_mod
from importlib import reload
reload(SIO_mod)

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
ENSO_data = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_data.xlsx')
ENSO_data_recent = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_recent_data.xlsx')
precip_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_Precip_data.csv')
PDO_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_PDO_data.csv', skiprows = 1)
path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

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

# create dataframe with spectra for each variable
spectra_temp_df = pd.DataFrame(columns = ['Temp_freq', 'Temp_spec', 'Temp_fft'])
spectra_sal_df = pd.DataFrame(columns = ['Sal_freq', 'Sal_spec', 'Sal_fft'])
spectra_PDO_df = pd.DataFrame(columns = ['PDO_freq', 'PDO_spec', 'PDO_fft'])
spectra_ENSO_df = pd.DataFrame(columns = ['ENSO_freq', 'ENSO_spec', 'ENSO_fft'])
spectra_tb_df = pd.DataFrame(columns = ['tb_freq', 'tb_spec', 'tb_fft'])

temp_sampled = np.mean(temp_output[0:37530].reshape(-1, 30), axis=1) #length = 1251
# compute spectral variables for each variable
for j in range(0,5):
	data_sets = [temp_data['SURF_TEMP_C_DETREND'], sal_data['SURF_SAL_PSU_DETREND'], PDO_data['Value'][744:], ENSO_data_all['VALUE'], temp_sampled]
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
	if j == 4:
		spectra_tb_df['tb_freq'] = freq
		spectra_tb_df['tb_spec'] = spec
		spectra_tb_df['tb_fft'] = fft
# find confidence intervals
spec_set = [spectra_temp_df['Temp_spec'], spectra_sal_df['Sal_spec'], spectra_PDO_df['PDO_spec'], spectra_ENSO_df['ENSO_spec'], spectra_tb_df['tb_spec']]
freq_set = [spectra_temp_df['Temp_freq'], spectra_sal_df['Sal_freq'], spectra_PDO_df['PDO_freq'], spectra_ENSO_df['ENSO_freq'], spectra_tb_df['tb_freq']]

df = 2
conf_lim = 0.95
conf_temp_l, conf_temp_h = SIO_mod.conf_int(spec_set[0], conf_lim, df)
conf_sal_l, conf_sal_h = SIO_mod.conf_int(spec_set[1], conf_lim, df)
conf_PDO_l, conf_PDO_h = SIO_mod.conf_int(spec_set[2], conf_lim, df)
conf_ENSO_l, conf_ENSO_h = SIO_mod.conf_int(spec_set[3], conf_lim, df) 
conf_tb_l, conf_tb_h = SIO_mod.conf_int(spec_set[4], conf_lim, 60) 

conf_l_set = [conf_temp_l, conf_sal_l, conf_PDO_l, conf_ENSO_l, conf_tb_l]
conf_h_set = [conf_temp_h, conf_sal_h, conf_PDO_h, conf_ENSO_h, conf_tb_h]

# color specifications
t_color = 'cadetblue'
s_color = 'darkslateblue'
p_color = 'seagreen'
e_color = 'steelblue'

# dictionaries to define legend names, colors, and labels
legend_name = ['Temp', 'Sal', 'PDO', 'ENSO', 'Band Averaged Temp']
color_sets = [t_color, s_color, p_color, e_color, t_color]
spectra_label = ['($^\circ$C$^2$)/(cycles/day)', '(PSU$^2$)/(cycles/day)','(PDO Index$^2$)/(cycles/day)','(ENSO Index$^2$)/(cycles/day)', '($^\circ$C$^2$)/(cycles/day)']
title_label = ['Power Spectrum SIO Temperature', 'Power Spectrum SIO Salinity', 'Power Spectrum PDO Index','Power Spectrum ENSO Index', 'Power Spectrum SIO Temeprature, Band Av = 30']
x_lab = ['Day', 'Day', 'Day', 'Day', 'Day']
freq_ann = 2*np.pi/365.25

# plot spectra for each variable with confidence intervals
for j in range(0,5):
	freq, spec, spec_amp, fft, delt, freq_T, freq_nyquist = SIO_mod.var_fft(data_sets[j])
	fig = plt.figure(figsize = (10,6))
	ax = fig.add_subplot(111)
	ax.loglog(freq_set[j], spec_set[j], color = color_sets[j])
	ax.fill_between(freq_set[j], conf_l_set[j], conf_h_set[j] , color = 'k', label = 'Conf. Lim. = ' + str(conf_lim), alpha = 0.2)
	ax.set_title(title_label[j])
	ax.set_ylabel('Energy Density \n' + spectra_label[j])
	ax.axvline(freq_nyquist, color = 'black', linestyle = '--', alpha = 0.5)
	ax.text(0.9, 0.9,'$\omega_{max}$', alpha = 0.5, transform = ax.transAxes)
	ax.axvline(freq_T, color = 'black', linestyle = '--', alpha = 0.5)
	ax.text(0.05, 0.9,'$\omega_o$', alpha = 0.5, transform = ax.transAxes)
	ax.axvline(freq_ann, color = 'black', linestyle = '--', alpha = 0.5)
	ax.legend(loc = 'lower left')
	if j <=1:
		ax.text(0.5, 0.9, 'Annual', alpha = 0.5, transform = ax.transAxes)
	if j >1:
		ax.text(0.72, 0.9, 'Annual', alpha = 0.5, transform = ax.transAxes)
	ax.set_xlabel('Cycles per ' + x_lab[j])
	if j <= 1:
		ax.set_ylim(10**-8, 10**5)
	im_name = title_label[j] + '.jpg'
	plt.savefig(path_out + im_name)
	plt.show()


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

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
ENSO_data = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_data.xlsx')
ENSO_data_recent = pd.read_excel('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_ENSO_recent_data.xlsx')
precip_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_Precip_data.csv')
path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

# convert year, month, day columns to single DATE column
sal_data['DATE'] = pd.to_datetime(sal_data[['YEAR', 'MONTH', 'DAY']])
temp_data['DATE'] = pd.to_datetime(temp_data[['YEAR', 'MONTH', 'DAY']])
precip_data['Date'] = pd.to_datetime(precip_data['DATE'], format='%Y-%m-%d')
ENSO_data_all = ENSO_data.append(ENSO_data_recent[323:], ignore_index = True)
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
temp_tri = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 60, min_periods = 3, win_type = 'triang').mean()

# plot temperature with ENSO
tstr = 'Temperature Anomaly and ENSO Index'
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(temp_data['DATE'], temp_data['SURF_TEMP_C_DETREND'], color = 'black', alpha = 0.5)
ax.set_ylabel('Temperature Anomaly ($^\circ$C)')
ax.set_xlabel('Date')
ax2 = ax.twinx()
ax2.fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']>0), color = 'r')
ax2.fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']<0), color = 'b')
ax2.set_ylabel('ENSO Index')
ax.set_title(tstr)
plt.show()

# butterworth low pass filter for temperature and salinity
fs = 1 # sampling frequency, once per day
fc = 1/500 # cut-off frequency of the filter (cut off periods shorter than 500 days)
w = fc / (fs / 2) #normalize the frequency
b, a = signal.butter(4, w, 'low')
temp_output = signal.filtfilt(b, a, temp_tri)
sal_output = signal.filtfilt(b, a, sal_tri)

t_color = 'cadetblue'
s_color = 'darkslateblue'

NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
# subplot 1, temperature and ENSO
axes[0].set_xlabel('Date')
axes[0].set_ylabel('Temperature Anomaly ($^\circ$C)')
axes[0].set_title('Temperature Anomalies, Butterworth filter cutoff = 500 days')
axes[0].plot(temp_data['DATE'],temp_output, color = t_color, linewidth = 2, label='Temp')
axes[0].set_ylabel('Temperature Anomaly ($^\circ$C)')
# plt.legend(loc = 'lower left')
# subplot 2, salinity and ENSO
axes[1].set_xlabel('Date')
axes[1].set_title('Salinity Anomalies, Butterworth filter cutoff = 500 days')
axes[1].plot(sal_data['DATE'],sal_output, color = s_color, linewidth = 2, label='Salinity')
axes[1].set_ylabel('Salinity Anomaly (PSU)')
# plt.legend(loc = 'lower left')
fig.tight_layout(pad=2.0)
im_name = 'TempSalAnomalies_NoSsn.jpg'
plt.savefig(path_out + im_name)
plt.show()

NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
# subplot 1, temperature and ENSO
axes[0].fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']>0), color = 'r', alpha = 0.6, label = '+ ENSO')
axes[0].fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']<0), color = 'b', alpha = 0.6, label = '- ENSO')
axes[0].set_xlabel('Date')
axes[0].set_ylabel('ENSO Index')
axes[0].set_title('Temperature Anomalies and ENSO Index')
ax2 = axes[0].twinx()
ax2.plot(temp_data['DATE'],temp_output, color = 'black', linewidth = 2, label='Temp')
ax2.set_ylabel('Temperature Anomaly ($^\circ$C)')
# axes[0].legend(loc = 'lower left')
# ax2.legend(loc = 'lower left')
plt.legend(loc = 'lower left')
# subplot 2, salinity and ENSO
axes[1].fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']>0), color = 'r', alpha = 0.6, label='+ ENSO')
axes[1].fill_between(ENSO_data_all['DATE'], ENSO_data_all['VALUE'], where=(ENSO_data_all['VALUE']<0), color = 'b', alpha = 0.6, label = '- ENSO')
axes[1].set_ylabel('ENSO Index')
axes[1].set_xlabel('Date')
axes[1].set_title('Salinity Anomalies and ENSO Index')
ax2 = axes[1].twinx()
ax2.plot(sal_data['DATE'],sal_output, color = 'black', linewidth = 2, label='Salinity')
ax2.set_ylabel('Salinity Anomaly (PSU)')
plt.legend(loc = 'lower left')
fig.tight_layout(pad=2.0)
im_name = 'ENSO_TempSalAnomalies.jpg'
plt.savefig(path_out + im_name)
plt.show()

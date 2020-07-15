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

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
PDO_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/NOAA_PDO_data.csv', skiprows = 1)
path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

# convert year, month, day columns to single DATE column
sal_data['DATE'] = pd.to_datetime(sal_data[['YEAR', 'MONTH', 'DAY']])
temp_data['DATE'] = pd.to_datetime(temp_data[['YEAR', 'MONTH', 'DAY']])
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
temp_tri = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 60, min_periods = 3, win_type = 'triang').mean()
PDO_ma = PDO_data['Value'][744:].rolling(center = True, window = 5, min_periods = 3, win_type = 'boxcar').mean()
PDO_ma2 = PDO_data['Value'][744:].rolling(center = True, window = 60, min_periods = 3).mean()


fig = plt.figure()
ax = fig.add_subplot(111)
# ax.plot(PDO_data['Date'][744:], PDO_data['Value'][744:], color = 'black', alpha = 0.5)
# ax.plot(PDO_data['Date'][744:], PDO_ma, color = 'black')
ax.plot(PDO_data['Date'][744:], PDO_ma2, color = 'red')
plt.show()

# remove seasonal temp
x = np.linspace(0,37538,37538)
y = 5*np.sin((2*np.pi/365)*x+(0.4005*365))
temp_data['SURF_TEMP_C_NOSSN'] = temp_tri - y

# plot temperature with PDO
# tstr = 'Temperature Anomaly and PDO Index'
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.plot(temp_data['DATE'], temp_data['SURF_TEMP_C_DETREND'], color = 'black', alpha = 0.5)
# ax.set_ylabel('Temperature Anomaly ($^\circ$C)')
# ax.set_xlabel('Date')
# ax2 = ax.twinx()
# ax2.fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]>0), color = 'r')
# ax2.fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]<0), color = 'b')
# ax2.set_ylabel('PDO Index')
# ax.set_title(tstr)
# plt.show()

# plot temperature anomaly and seasonal sinusoid fit; plot the difference
# NR = 2; NC = 1
# fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
# axes[0].plot(temp_data['DATE'], temp_ma, color = 'magenta', label = 'Smoothed Data')
# axes[0].plot(temp_data['DATE'], y, 'orange', label = 'Seasonal Sinusoid Fit')
# axes[0].set_xlabel('Date')
# axes[0].set_ylabel('Temperature Anomaly ($^\circ$C)')
# axes[0].set_title('Temperature Anomaly and Seasonal Fit')
# axes[0].legend(loc = 'upper left')
# axes[1].plot(temp_data['DATE'], temp_data['SURF_TEMP_C_NOSSN'], color = 'green')
# axes[1].set_xlabel('Date')
# axes[1].set_ylabel('Temperature Anomaly ($^\circ$C)')
# axes[1].set_title('Temperature Anomaly, Seasons Removed')
# fig.tight_layout(pad=2.0)
# plt.show()

# butterworth low pass filter for temperature
fs = 1 # sampling frequency, once per day
fc = 1/500 # cut-off frequency of the filter
w = fc / (fs / 2) #normalize the frequency
b, a = signal.butter(4, w, 'low')
temp_output = signal.filtfilt(b, a, temp_data['SURF_TEMP_C_DETREND'])
sal_output = signal.filtfilt(b, a, sal_data['SURF_SAL_PSU_DETREND'])

NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
# subplot 1, temperature and PDO
axes[0].fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]>0), color = 'r', alpha = 0.6, label = '+ PDO')
axes[0].fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]<0), color = 'b', alpha = 0.6, label = '- PDO')
axes[0].set_xlabel('Date')
axes[0].set_ylabel('PDO Index')
axes[0].set_title('Temperature Anomalies and PDO Index')
ax2 = axes[0].twinx()
ax2.plot(temp_data['DATE'],temp_output, color = 'black', linewidth = 2, label='Temp')
ax2.set_ylabel('Temperature Anomaly ($^\circ$C)')
# axes[0].legend(loc = 'lower left')
# ax2.legend(loc = 'lower left')
plt.legend(loc = 'lower left')
# subplot 2, precipitation and PDO
axes[1].fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]>0), color = 'r', alpha = 0.6, label='+ PDO')
axes[1].fill_between(PDO_data['DATE'][744:], PDO_data['Value'][744:], where=(PDO_data['Value'][744:]<0), color = 'b', alpha = 0.6, label = '- PDO')
axes[1].set_ylabel('PDO Index')
axes[1].set_xlabel('Date')
axes[1].set_title('Salinity Anomalies and PDO Index')
ax2 = axes[1].twinx()
ax2.plot(sal_data['DATE'],sal_output, color = 'black', linewidth = 2, label='Salinity')
ax2.set_ylabel('Salinity Anomaly (PSU)')
plt.legend(loc = 'lower left')
fig.tight_layout(pad=2.0)
im_name = 'PDO_TempSalAnomalies.jpg'
plt.savefig(path_out + im_name)

plt.show()

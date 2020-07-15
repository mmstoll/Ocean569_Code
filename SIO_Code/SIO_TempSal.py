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
import SIO_modules as SIO_mod
from importlib import reload
reload(SIO_mod)

# read in temp and sal files
sal_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_SALT_1916-201905.txt', sep='\t', skiprows = 27)
temp_data = pd.read_csv('/Users/MMStoll/Python/Data/Ocean569_Data/SIO_Data/SIO_TEMP_1916_201905.txt', sep='\t', skiprows = 26)
path_out = '/Users/MMStoll/Python/Output/Ocean569_Output/SIO_Output/'

# convert year, month, day columns to single DATE column
sal_data['DATE'] = pd.to_datetime(sal_data[['YEAR', 'MONTH', 'DAY']])
temp_data['DATE'] = pd.to_datetime(temp_data[['YEAR', 'MONTH', 'DAY']])

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

sal_tri = sal_data['SURF_SAL_PSU_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'triang').mean()
temp_tri = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 60, min_periods = 3, win_type = 'triang').mean()

t_color = 'cadetblue'
s_color = 'darkslateblue'
p_color = 'seagreen'
e_color = 'steelblue'

# plot the time series, no average
tstr = 'SIO Temperature and Salinity Time Series'
im_name = 'SIO_TempSal_TimeSeries.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].plot(temp_data['DATE'],temp_data['SURF_TEMP_C_NOAVG'], color=t_color, alpha = 0.5)
axes[0].plot(temp_data['DATE'],temp_fit_fn(temp_data.index), color = t_color, linewidth = 3)
axes[0].set_xlabel('Date')
axes[0].set_ylabel('Temperature ($^\circ$C)')
axes[1].plot(sal_data['DATE'],sal_data['SURF_SAL_PSU_NOAVG'], color=s_color, alpha = 0.5)
axes[1].plot(sal_data['DATE'],sal_fit_fn(sal_data.index), color =s_color, linewidth = 3)
axes[1].set_xlabel('Date')
axes[1].set_ylabel('Salinity (PSU)')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

# plot the time series with the trend and average removed
tstr = 'SIO Temperature and Salinity Time Series Anomalies'
im_name = 'SIO_TempSal_TimeSeries_Anomalies.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].plot(temp_data['DATE'],temp_data['SURF_TEMP_C_DETREND'], color=t_color)
axes[0].set_xlabel('Date')
axes[0].set_ylabel('Temperature Anomaly ($^\circ$C)')
axes[1].plot(sal_data['DATE'],sal_data['SURF_SAL_PSU_DETREND'], color=s_color)
axes[1].set_xlabel('Date')
axes[1].set_ylabel('Salinity Anomaly (PSU)')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

# create variables to carry out spectral analysis
ss = len(sal_data['SURF_SAL_PSU_NOAVG'])
tt = len(temp_data['SURF_TEMP_C_NOAVG'])
pi = np.pi
ss_half = int(ss/2+1)
tt_half = int(tt/2+1)
ss_freq = np.zeros(ss_half)
tt_freq = np.zeros(tt_half)

ss_T_length = ss
ss_delt = 1 #days
tt_T_length = tt
tt_delt = 1 #days

ss_freq_nyquist=pi/ss_delt
ss_freq_T=2.*pi/(ss*ss_delt)
tt_freq_nyquist=pi/tt_delt
tt_freq_T=2.*pi/(tt*tt_delt)  
freq_ann = 2*pi/365

ss_omega0 = 2.*pi/(ss_T_length*ss_delt)
tt_omega0 = 2.*pi/(tt_T_length*tt_delt)
for i in range(0,ss_half):
    ss_freq[i]=i*ss_omega0
for i in range(0,tt_half):
	tt_freq[i] = i*tt_omega0

# compute the power spectrum for temperature and salinity, detrended
ss_fft = np.fft.rfft(sal_data['SURF_SAL_PSU_DETREND'],n=ss)
ss_fourier_amp = np.sqrt((np.real(ss_fft)**2+np.imag(ss_fft)**2))
ss_fourier_phase=180.*np.arctan2(np.imag(ss_fft),np.real(ss_fft))/pi
ss_spec = np.real(ss_fft*np.conj(ss_fft))/(2.*pi*ss_T_length*ss_delt)
ss_spec_amp = (np.absolute(ss_fft))**2/(2*pi*ss_T_length*ss_delt)

tt_fft = np.fft.rfft(temp_data['SURF_TEMP_C_DETREND'],n=tt)
tt_fourier_amp = np.sqrt((np.real(tt_fft)**2+np.imag(tt_fft)**2))
tt_fourier_phase=180.*np.arctan2(np.imag(tt_fft),np.real(tt_fft))/pi
tt_spec = np.real(tt_fft*np.conj(tt_fft))/(2.*pi*tt_T_length*tt_delt)
tt_spec_amp = (np.absolute(tt_fft))**2/(2*pi*tt_T_length*tt_delt)

# plot the spectra, detrended
tstr = 'SIO Temperature and Salinity Power Spectra'
im_name = 'SIO_TempSal_Spectra.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].loglog(tt_freq,tt_spec_amp, color=t_color)
axes[0].plot([tt_freq_nyquist,tt_freq_nyquist], [10**-8,10**6],'--k', alpha = 0.5)
axes[0].text(10**0.16,10**4,'$\omega_{max}$', alpha = 0.5)
axes[0].plot([tt_freq_T,tt_freq_T], [10**-8,10**6], '--k', alpha = 0.5)
axes[0].text(10**-3.7,10**4,'$\omega_o$', alpha = 0.5)
axes[0].plot([freq_ann,freq_ann], [10**-8,10**6],'--k', alpha = 0.5)
axes[0].set_ylim(10**-7, 10**5)
axes[0].set_xlabel('$\omega$ (radians/day)')
axes[0].set_ylabel('Energy Density \n($^\circ$C$^2$)/(radians/day)')
axes[1].loglog(ss_freq,ss_spec_amp, color=s_color)
axes[1].plot([ss_freq_nyquist,ss_freq_nyquist], [10**-10,10**2], '--k', alpha = 0.5)
axes[1].text(10**0.16,10**1,'$\omega_{max}$', alpha = 0.5)
axes[1].plot([ss_freq_T,ss_freq_T], [10**-10,10**2], '--k', alpha = 0.5)
axes[1].text(10**-3.7,10**1,'$\omega_o$', alpha = 0.5)
axes[1].plot([freq_ann,freq_ann], [10**-10,10**2],'--k', alpha = 0.5)
axes[1].set_ylim(10**-9, 10**2)
axes[1].set_xlabel('$\omega$ (radians/day)')
axes[1].set_ylabel('Energy Density \n(PSU$^2$)/(radians/day)')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

# filter data
sal_ma = sal_data['SURF_SAL_PSU_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'boxcar').mean()
sal_tri = sal_data['SURF_SAL_PSU_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'triang').mean()
temp_ma = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'boxcar').mean()
temp_tri = temp_data['SURF_TEMP_C_DETREND'].rolling(center = True, window = 30, min_periods = 3, win_type = 'triang').mean()

# plot the temperature time series, running mean and triangle filters
tstr = 'SIO Temperature Boxcar and Triangle Filters'
im_name = 'SIO_TempTime_Filters.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].plot(temp_data['DATE'],temp_data['SURF_TEMP_C_DETREND'], color = t_color, alpha = 0.5)
axes[0].plot(temp_data['DATE'],temp_ma, color=t_color, linewidth = 1, label = 'Boxcar Interval = 30 days')
axes[0].set_xlabel('Date')
axes[0].set_ylabel('Temperature ($^\circ$C)')
axes[0].legend(loc='lower left')
axes[1].plot(temp_data['DATE'],temp_data['SURF_TEMP_C_DETREND'], color = t_color, alpha = 0.5)
axes[1].plot(temp_data['DATE'],temp_tri, color=t_color, linewidth = 1, label = 'Triangle Interval = 30 days')
axes[1].set_xlabel('Date')
axes[1].set_ylabel('Temperature ($^\circ$C)')
axes[1].legend(loc='lower left')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

# plot the salinity time series, running mean and triangle filters
tstr = 'Scripps Pier Triangle Filtered Data'
im_name = 'SIO_SalTime_Filters.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].plot(temp_data['DATE'],temp_data['SURF_TEMP_C_DETREND'], color = t_color, alpha = 0.5)
axes[0].plot(temp_data['DATE'],temp_tri, color=t_color, linewidth = 1, label = 'Triangle Interval = 10 days')
axes[0].set_xlabel('Date')
axes[0].set_ylabel('Temperature ($^\circ$C)')
axes[0].legend(loc='lower left')
axes[1].plot(sal_data['DATE'],sal_data['SURF_SAL_PSU_DETREND'], color = s_color, alpha = 0.5)
axes[1].plot(sal_data['DATE'],sal_tri, color=s_color, linewidth = 1, label = 'Triangle Interval = 10 days')
axes[1].set_xlabel('Date')
axes[1].set_ylabel('Salinity (PSU)')
axes[1].legend(loc='lower left')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

def band_average(fft_var1,fft_var2,frequency,n_av):
	# fft_var1 and fft_var2 are the inputs computed via fft
	# they can be the same variable or different variables
	# n_av is the number of bands to be used for smoothing (nice if it is an odd number)
	# this function is limnited to 100,000 points but can easily be modified
    nmax=100000

	# define some variables and arrays
    n_spec=len(fft_var1)
    n_av2=int(n_av//2+1)
    spec_amp_av=np.zeros(nmax)
    spec_phase_av=np.zeros(nmax)
    freq_av=np.zeros(nmax)
	# average the lowest frequency bands first (with half as many points in the average)
    sum_low_amp=0.
    sum_low_phase=0.
    count=0
    spectrum_amp=np.absolute(fft_var1*np.conj(fft_var2))
    spectrum_phase=np.angle(fft_var1*np.conj(fft_var2),deg=True)
	#
    for i in range(0,n_av2):
        sum_low_amp+=spectrum_amp[i]
        sum_low_phase+=spectrum_phase[i]
    spec_amp_av[0]=sum_low_amp/n_av2
    spec_phase_av[0]=sum_low_phase/n_av
	# compute the rest of the averages
    for i in range(n_av2,n_spec-n_av,n_av):
        count+=1
        spec_amp_est=np.mean(spectrum_amp[i:i+n_av])
        spec_phase_est=np.mean(spectrum_phase[i:i+n_av])
        freq_est=frequency[i+n_av//2]
        spec_amp_av[count]=spec_amp_est
        spec_phase_av[count]=spec_phase_est
        freq_av[count]=freq_est
	# contract the arrays
    spec_amp_av=spec_amp_av[0:count]
    spec_phase_av=spec_phase_av[0:count]
    freq_av=freq_av[0:count]
    return spec_amp_av,spec_phase_av,freq_av,count

# band average the salinity, temperature, and coherence 
n_av = 30
tt_fft_star = np.conj(tt_fft)
sal_spec,sal_phase,ss_freq_av,count=band_average(ss_fft,ss_fft,ss_freq,n_av)
temp_spec,temp_phase,temp_freq_av,count=band_average(tt_fft,tt_fft,tt_freq,n_av)
cospec_amp,cospec_phase,freq_av,count=band_average(ss_fft,tt_fft_star,tt_freq,n_av)
coh_sq=cospec_amp**2/(temp_spec*sal_spec)

# plot the coherence and phase between salinity and temperature
tstr = 'SIO Temperature and Salinity Coherence and Phase'
im_name = 'SIO_TempSal_CoherencePhase.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].semilogx(freq_av,coh_sq, color = 'black')
axes[0].set_xlabel('$\omega$ (radians/day)')
axes[0].set_ylabel('Squared Coherence $\it{T}$-$\it{S}$')
# axes[0].grid()
axes[1].semilogx(freq_av, cospec_phase, color = 'black')
axes[1].set_xlabel('$\omega$ (radians/day)')
axes[1].set_ylabel('Phase $\it{T}$-$\it{S}$, degrees')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()

n_av = 30
# plot the coherence and phase between salinity and temperature after smoothing 
t_freq,t_spec,t_spec_amp,t_fft,t_delt,t_freq_T,t_freq_nyquist = SIO_mod.var_fft(temp_tri)
s_freq,s_spec,s_spec_amp,s_fft,s_delt,s_freq_T,s_freq_nyquist = SIO_mod.var_fft(sal_tri)
s_spec_b,s_phase_b,s_freq_av_b,count=band_average(s_fft,s_fft,s_freq,n_av)
t_spec_b,t_phase_b,t_freq_av_b,count=band_average(t_fft,t_fft,t_freq,n_av)
s_fft_star = np.conj(s_fft)
cospec_amp2,cospec_phase2,freq_av2,count2=band_average(t_fft,s_fft_star,t_freq,n_av)
coh_sq2=cospec_amp2**2/(t_spec_b*s_spec_b)

# plot the coherence and phase between salinity and temperature
tstr = 'SIO Temperature and Salinity Coherence and Phase'
im_name = 'SIO_TempSal_CoherencePhase_små.jpg'
NR = 2; NC = 1
fig, axes = plt.subplots(nrows = NR,ncols=NC,figsize = (10,6))
axes[0].semilogx(freq_av2,coh_sq2, color = 'r')
axes[0].set_xlabel('$\omega$ (radians/day)')
axes[0].set_ylabel('Squared Coherence $\it{T}$-$\it{S}$')
# axes[0].grid()
axes[1].semilogx(freq_av2, cospec_phase2, color = 'r')
axes[1].set_xlabel('$\omega$ (radians/day)')
axes[1].set_ylabel('Phase $\it{T}$-$\it{S}$, degrees')
fig.suptitle(tstr)
plt.savefig(path_out + im_name)
plt.show()



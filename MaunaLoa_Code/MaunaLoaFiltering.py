#program to create spectra of MaunaLoa CO2 data
import numpy as np
import matplotlib.pyplot as plt
import random
from netCDF4 import Dataset
#
# read and plot Mauna Loa CO2 data, from a netcdf file
# 
# define the netcdf reading function
def import_data(file_name):
    data_netcdf = Dataset(file_name, mode = 'r')
    data = {}
    for vname in list(data_netcdf.variables):
        data[str(vname)] = data_netcdf.variables[vname][:]
    data_netcdf.close()
    return data
#
path_in = '/Users/MMStoll/Python/Data/Ocean569_Data/MaunaLoa_Data/MaunaLoa.nc'
ML_data = import_data(path_in)
#
nn=len(ML_data['T'])
mm=int(nn/2+1)
#
# define the data arrays
time_meas=np.zeros(nn)
co2_meas=np.zeros(nn)
cc=np.zeros(nn)
cc_detrend=np.zeros(nn)
freq=np.zeros(mm)
r=np.zeros(nn)
#
#parse the time, CO2 data from the input file
for i in range(0,nn):
    time_meas[i]=(float(ML_data['T'][i])-float(ML_data['T'][0]))
    co2_meas[i]=float(ML_data['co2'][i])
#
## find the value of CO2 mean 
co2mean=np.nanmean(co2_meas) #remove nans from mean!!
print('The co2 mean is:')
print(co2mean)
co2_meas_nonan = co2_meas[~np.isnan(co2_meas)]
# 
# remove the CO2 mean from the data
for i in range(0,nn):
  cc[i]=co2_meas[i]-co2mean
# 
# remove NaN values from the dataset  
cc_nonan = cc[~np.isnan(cc)] #remove NaN from dataset!! 
# 
#detrend the CO2 data, remove increase over time 
for i in range(0,nn):
  r[i] = .1066*i+310.66
  cc_detrend[i]=co2_meas[i]-r[i] #remove line from CO2_meas data 
cc_nonan_detrend = cc_detrend[~np.isnan(cc_detrend)]
# determine the frequencies for the spectrum
delt=1
T_length=nn
pi=np.pi
omega0=2.*pi/(T_length*delt)
#
for i in range(0,mm):
     freq[i]=i*omega0
#         
# compute the fft of the input data (temperature in this case)
# the fft will yield a set of complex numbers
# also compute their complex conjugates
# multiply these together to get the spectrum
zz=np.fft.rfft(cc_nonan_detrend,n=nn)
fourier_amp_detrend=np.sqrt((np.real(zz)**2+np.imag(zz)**2))
fourier_phase=180.*np.arctan2(np.imag(zz),np.real(zz))/pi
spec_cc_detrend=np.real(zz*np.conj(zz))/(2.*pi*T_length*delt)
spec_cc_amp_detrend=(np.absolute(zz))**2/(2*pi*T_length*delt)

# =================================================================================
n_len = 5
omega_chk = 2.*pi*nn*(1./n_len)
amp = 1./(n_len + 1)

tt = np.zeros(nn)
yy = np.zeros(nn)
freq0=np.zeros(mm)
freq1=np.zeros(mm)
sum_signal = 0.

# ================================================================================
# Running Mean
for i in range(0,nn):
  tt[i] = i
  if i <= n_len//2:
    signal = amp
  else:
    signal = 0.
  yy[i]=signal
  sum_signal+=yy[i]

for i in range(0,mm):
  freq0[i]=i*omega0
  freq1[i]=freq0[i]*T_length

# carry out the fft
ww = np.fft.rfft(yy,n=nn)
ww_real=ww.real
ww_imag=ww.imag
ww_mag=np.sqrt(ww_real**2+ww_imag**2)

fig1=plt.figure()
plt.xlabel('Frequency')
plt.title('Running Mean Filter (5 months)')
# plt.ylabel()
# plt.semilogx(freq1,ww_mag)
# filtered_running_mean = ww_mag*fourier_amp_detrend
filtered_running_mean = ww_mag*fourier_amp_detrend
plt.loglog(freq,filtered_running_mean, color='blue', label='FT Filtered Data')
plt.loglog(freq,fourier_amp_detrend, color='red', label='FT Actual Data')
plt.xlabel('$\omega$ (radians/month)',fontsize=15,ha='center')
plt.ylabel('Fourier Amplitude (ppm)',fontsize=15)
plt.legend(loc=3)
plt.show()

fig2 = plt.figure()
zz_detrend = np.fft.irfft(filtered_running_mean,n=nn)
# zz_trend=np.fft.rfft(co2_meas_nonan,n=nn)
# fourier_amp_trend=np.sqrt((np.real(zz_trend)**2+np.imag(zz_trend)**2))
# filtered_running_mean_trend = ww_mag*fourier_amp_trend
# filtered_timeseries = np.fft.irfft(filtered_running_mean_trend,n=nn)
plt.title('Time Series - Filtered and Actual Data')
plt.plot(time_meas,zz_detrend, color ='blue', label='Filtered Data')
plt.plot(time_meas,cc_detrend, color='red', label='Actual Data')
plt.ylabel('CO2 level Anomalies (ppm)')
plt.legend(loc=3)
plt.xlabel('Time (Months since 1956)')
plt.show()
# ================================================================================
# Triangle Filter
n_len_tri = 5
omega_chk_tri = 2.*pi*nn*(1./n_len_tri)
amp_tri = 1

tt = np.zeros(nn)
yy_tri = np.zeros(nn)
freq0=np.zeros(mm)
freq1=np.zeros(mm)
sum_signal_tri = 0.

for i in range(0,nn):
  if i <= n_len_tri//2:
    signal_tri=i*4./(n_len_tri)**2
  if i > n_len_tri//2 and i <= n_len_tri:
    signal_tri=(2./n_len_tri)*(1.-(i-n_len_tri/2)*4/((n_len_tri)**2))
  if i > n_len_tri:
    signal_tri=0.
  yy_tri[i]=signal_tri
  sum_signal+=yy_tri[i]

# carry out the fft
uu = np.fft.rfft(yy,n=nn)
uu_real=uu.real
uu_imag=uu.imag
uu_mag=np.sqrt(uu_real**2+uu_imag**2)

fig3=plt.figure()
filtered_tri = uu_mag*fourier_amp_detrend
plt.title('Triangle Filter (5 months)')
plt.xlabel('$\omega$ (radians/month)',fontsize=15,ha='center')
plt.ylabel('Fourier Amplitude (ppm)',fontsize=15)
plt.loglog(freq,filtered_tri, color='blue', label='Filtered FT Data')
plt.loglog(freq,fourier_amp_detrend, color='red', label='FT Data')
plt.legend(loc=3)
plt.show()

# function in pandas module that does these filters for you in the time domain

#OOI HOT Mooring Data 
#file to read and plot Fourier amplitudes of temperature
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import math
#
# read and plot Mauna Loa CO2 data, from a netcdf file
#
# define the netcdf reading function
def import_data(file_name):   
  #'''This function works regardless of the file called as long as a list of variables can be identified easily in the file.''' 
  data_netcdf = Dataset(file_name, mode = 'r')
  data = {}
    # 
  for vname in list(data_netcdf.variables):
    data[str(vname)] = data_netcdf.variables[vname][:]
  data_netcdf.close()
  return data
#
# main program
#
# define the input and output files
path_in='/Users/MMStoll/Python/Data/Ocean569_Data/MaunaLoa_Data/MaunaLoa.nc'
path_out1='/Users/MMStoll/Python/Output/Ocean569_Output/MaunaLoa_Output/MaunaLoaPlots_Output/MaunaLoa_time.plot.jpg'
path_out2='/Users/MMStoll/Python/Output/Ocean569_Output/MaunaLoa_Output/MaunaLoaPlots_Output/MaunaLoa_time_nomean.plot.jpg'
path_out3='/Users/MMStoll/Python/Output/Ocean569_Output/MaunaLoa_Output/MaunaLoaPlots_Output/MaunaLoa_time_detrend.plot.jpg'
path_out4='/Users/MMStoll/Python/Output/Ocean569_Output/MaunaLoa_Output/MaunaLoaPlots_Output/MaunaLoa_freq_nomean.plot.jpg'
path_out5='/Users/MMStoll/Python/Output/Ocean569_Output/MaunaLoa_Output/MaunaLoaPlots_Output/MaunaLoa_freq_detrend.plot.jpg'
#path_out2='/Users/MMStoll/Documents/Python/Ocean569/HOTexample/HOT_fourier.amplitudes.jpg'
#path_out3='/Users/MMStoll/Documents/Python/Ocean569/HOTexample/full_averaged_spectra.jpg'

# read the input file (netcdf)
ML_data=import_data(path_in)

# determine the length of the data and the resulting number of Fourier estimates
nn=len(ML_data['T'])
mm=int(nn/2+1)
oo=len(ML_data['co2'])
pi=np.pi
#
# define the data arrays
time_meas=np.zeros(nn)
co2_meas=np.zeros(nn)
cc=np.zeros(nn)
cc_detrend=np.zeros(nn)
freq=np.zeros(mm)
r=np.zeros(nn)
sinfn=np.zeros(nn)
#
for i in range (0,nn):
  sinfn[i] = 3.5*math.sin(i*(2*pi)/12.)

#parse the time, CO2 data from the input file
for i in range(0,nn):
    time_meas[i]=(float(ML_data['T'][i])-float(ML_data['T'][0]))
    co2_meas[i]=float(ML_data['co2'][i])

#plot the CO2 time series
fig1=plt.figure()
#plt.xlim(0,400)
plt.plot(time_meas,co2_meas,color='firebrick')
plt.grid()
plt.xlabel('Months since Jan. 1960 (through 1995)')
plt.ylabel('CO2 Levels (ppm)')
plt.title('Mauna Loa Time Series')
plt.savefig(path_out1)
#plt.show()

# find the value of CO2 mean 
co2mean=np.nanmean(co2_meas) #remove nans from mean!!
print('The co2 mean is:')
print(co2mean)

# remove the CO2 mean from the data
for i in range(0,nn):
  cc[i]=co2_meas[i]-co2mean

# remove NaN values from the dataset  
cc_nonan = cc[~np.isnan(cc)] #remove NaN from dataset!!

#detrend the CO2 data, remove increase over time 
for i in range(0,nn):
  r[i] = .1066*i+310.66
  cc_detrend[i]=co2_meas[i]-r[i] #remove line from CO2_meas data
cc_nonan_detrend = cc_detrend[~np.isnan(cc_detrend)]

#plot time series with mean removed
fig2=plt.figure()
plt.plot(time_meas,cc)
plt.xlabel('Months since Jan. 1960 (through 1995)')
plt.ylabel('CO2 Levels (ppm from mean)')
plt.title('Mauna Loa Time Series (mean removed)')
plt.grid()
plt.show()
plt.savefig(path_out2)

#plot time series with trend removed
fig3=plt.figure()
plt.plot(time_meas,cc_detrend)
plt.xlabel('Months since Jan. 1960 (through 1995)')
plt.ylabel('CO2 Levels (ppm, detrended)')
plt.title('Mauna Loa Time Series (detrended)')
plt.grid()
plt.show()
plt.savefig(path_out3)
plt.plot(time_meas,sinfn,color='purple')

# plot time series with trend and mean removed 

# determine the frequencies for the spectrum
delt=1
T_length=nn
omega0=2.*pi/(T_length*delt)
# #
for i in range(0,mm):
     freq[i]=i*omega0
# #
# # compute the fft of the input data (temperature in this case)
# # the fft will yield a set of complex numbers

zz=np.fft.rfft(cc_nonan,n=nn)
fourier_amp_trend=np.sqrt((np.real(zz)**2+np.imag(zz)**2))
fourier_phase=180.*np.arctan2(np.imag(zz),np.real(zz))/pi

yy=np.fft.rfft(sinfn,n=nn)
yy_fourier_amp=np.sqrt((np.real(yy)**2+np.imag(yy)**2))
yy_fourier_phase=180.*np.arctan2(np.imag(yy),np.real(yy))/pi

# plot the Fourier amplitudes as a function of frequency
# show the maximum and minimum frequencies
fig4=plt.figure()
plt.ylim(10**-1,10**4)
plt.xlim(.0105665,4)
plt.loglog(freq,fourier_amp_trend,color='magenta')
plt.xlabel('$\omega$ (radians/month)',fontsize=15,ha='center')
plt.ylabel('Fourier Amplitude (ppm)',fontsize=15)
plt.title('Frequency Domain - Mauna Loa CO2 Levels')
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt) 
plt.plot([freq_nyquist,freq_nyquist],[0.01,3.e5],'--k')
plt.text(1.8,0.5,'$\omega_{max}$',fontsize=12,color='blue')
plt.plot([freq_T,freq_T],[0.01,3.e5],'--k',zorder=10)
plt.text(0.015,0.5,'$\omega_o$',fontsize=12,color='blue')
plt.grid()
plt.plot(freq,yy_fourier_amp,color='purple')
plt.savefig(path_out4)

# show the annual CO2 cycle
period_ann=12
freq_ann=2.*pi/period_ann
plt.text(freq_ann,732.8,'annual',color='blue')
#show the semi-annual CO2 cycle
period_semi_ann=6
freq_semi_ann=2.*pi/period_semi_ann
plt.text(freq_semi_ann,220,'semi-annual',color='blue')
# plt.savefig(path_out2)

#cycle every 4 months
period_4mon=3.9989
freq_4mon=2.*pi/period_4mon
#plt.text(freq_4mon,32,'4 months',color='blue')
#cycle every 3 months
period_3mon=2.9921
freq_3mon=2.*pi/period_3mon
#plt.text(freq_3mon,22,'3 months',color='blue')

#plot detrended Fourier Transform
zz_detrend=np.fft.rfft(cc_nonan_detrend,n=nn)
fourier_amp_detrend=np.sqrt((np.real(zz_detrend)**2+np.imag(zz_detrend)**2))
fourier_phase_detrend=180.*np.arctan2(np.imag(zz_detrend),np.real(zz_detrend))/pi

fig5=plt.figure()
plt.loglog(freq,fourier_amp_detrend,color='magenta')
plt.ylim(0.4,751)
plt.plot(freq,yy_fourier_amp,color='purple')

freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt) 
plt.plot([freq_nyquist,freq_nyquist],[0.4,750],'--k')
plt.text(1.8,0.5,'$\omega_{max}$',fontsize=12,color='blue')
plt.plot([freq_T,freq_T],[.4,750],'--k',zorder=10)
plt.text(0.015,0.5,'$\omega_o$',fontsize=12,color='blue')
plt.text(freq_ann,594,'annual',color='blue')
plt.text(freq_semi_ann,135,'semi-annual',color='blue')
plt.text(freq_4mon,29,'4 months',color='blue')
plt.text(freq_3mon,21,'3 months',color='blue')

plt.xlabel('$\omega$ (radians/month)',fontsize=15,ha='center')
plt.ylabel('Fourier Amplitude (ppm, detrended)',fontsize=15)
plt.title('Frequency Domain - Mauna Loa CO2 Levels detrended')
plt.grid()
plt.savefig(path_out5)
plt.show()

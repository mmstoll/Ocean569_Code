#OOI HOT Mooring Data 
file to read and plot Fourier amplitudes of temperature
#
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from netCDF4 import Dataset
#
# read and plot OOI HOT mooring data, from a netcdf file
#
# define the netcdf reading function
#
def import_data(file_name):   
    '''This function works regardless of the file called as long as a list of 
    variables can be identified easily in the file.'''
    
    data_netcdf = Dataset(file_name, mode = 'r')
    data = {}
#    
    for vname in list(data_netcdf.variables):
        data[str(vname)] = data_netcdf.variables[vname][:] 
#            
    data_netcdf.close()
#            
    return data
#
# main program
#
# define the input and output files
#
path_in='/Users/riser/Desktop/ocean.569A/OS_WHOTS_201606_D_MICROCAT-025m.nc'
path_out1='/Users/riser/Desktop/ocean.569A/HOT_data.plot.jpg'
path_out2='/Users/riser/Desktop/ocean.569A/HOT_fourier.amplitudes.jpg'
path_out3='/Users/riser/Desktop/ocean.569A/full_averaged_spectra.jpg'
#
# read the input file (netcdf)
#
HOT_data=import_data(path_in)
#
# determine the length of the data and the resulting number of Fourier estimates
#
nn=len(HOT_data['TIME'])
mm=int(nn/2+1)
#
# define the data arrays
#
time_meas=np.zeros(nn)
time_meas_days=np.zeros(nn)
temp_meas=np.zeros(nn)
sal_meas=np.zeros(nn)
tt=np.zeros(nn)
ss=np.zeros(nn)
freq=np.zeros(mm)
#
# parse the time, temperature, and salinity data from the input file

for i in range(0,nn):
   time_meas[i]=24.*(float(HOT_data['TIME'][i])-float(HOT_data['TIME'][0]))
   time_meas_days[i]=(float(HOT_data['TIME'][i])-float(HOT_data['TIME'][0]))
   temp_meas[i]=float(HOT_data['TEMP'][i])
   sal_meas[i]=float(HOT_data['PSAL'][i])
#
# plot the temperature time series
# 
fig1=plt.figure()
plt.xlim(0,400)
plt.plot(time_meas_days,temp_meas,color='firebrick')
plt.grid()
plt.xlabel('Days')
plt.ylabel('Temperature  ($^o$C)')
plt.title('HOT Microcat 2016 (3-min samples)')
#plt.savefig(path_out1)
plt.show()
#
# remove the temperature and salinity means from the data
#
tmean=np.mean(temp_meas)
smean=np.mean(sal_meas)
for i in range(0,nn):
    tt[i]=temp_meas[i]-tmean
    ss[i]=sal_meas[i]-smean
# 
# determine the frequencies for the spectrum
#
delt=0.00208333
T_length=nn
pi=np.pi
omega0=2.*pi/(T_length*delt)
#
for i in range(0,mm):
    freq[i]=i*omega0
#
# compute the fft of the input data (temperature in this case)
# the fft will yield a set of complex numbers
#
zz=np.fft.rfft(tt,n=nn)
fourier_amp=np.sqrt((np.real(zz)**2+np.imag(zz)**2))
fourier_phase=180.*np.arctan2(np.imag(zz),np.real(zz))/pi
#
# plot the Fourier amplitudes as a function of frequency
# show the maximum and minimum frequencies
#
fig2=plt.figure(figsize=(7,7))
plt.ylim(0.01,3.e5)
plt.loglog(freq,fourier_amp,color='magenta')
plt.xlabel('$\omega$ (radians/day)',fontsize=15,ha='center')
plt.ylabel('Fourier Amplitude ($^o$C)',fontsize=15)
freq_nyquist=pi/delt
freq_T=2.*pi/(nn*delt)
plt.plot([freq_nyquist,freq_nyquist],[0.01,3.e5],'--k')
plt.text(5.e2,5.e4,'$\omega_{max}$',fontsize=12,color='blue')
plt.plot([freq_T,freq_T],[0.01,3.e5],'--k',zorder=10)
plt.text(0.02,0.3,'$\omega_o$',fontsize=12,color='blue')
#
# show the K1 tide and add a grid
#
period_K1=23.93
freq_K1=2.*pi+0.8
plt.text(5.,2.5e3,'K1',color='blue')
plt.grid()
#plt.savefig(path_out2)
plt.show()
#
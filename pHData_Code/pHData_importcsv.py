import csv
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import datetime as dt
from sklearn.linear_model import LinearRegression

#read in csv file, skip the first row with headers, and rename the column titles with 'names'
#this command, pd.read_csv, automatically imports the data into a DataFrame
data = pd.read_csv('/Users/MMStoll/Documents/Python/Data/Ocean569/pHData/OOI_CE04OSSM-H1_pH.csv', skiprows=0, names=['DateTime', 'Depth', 'pH'])
#print(data.DateTime.dtype) #see what datatype DateTime is

#rearrage the order of the columns using 'columns' followed by preferred order
# data2 = pd.DataFrame(data, columns=['Depth', 'DateTime', 'pH'])
# print(data2)

#convert the columns, pH and DateTime, into a Series
pH = data.pH
date = data.DateTime

#unfortunately this doesnt work to convert a series into datetime format
# nn=len(date)
# date2=np.zeros(nn)
# date_format= "%Y-%m-%dT%H:%M:%S-%f"
# for i in range(0,nn):
#    date2[i]=datetime.strptime(date[i], date_format)

#convert time stamp series 'date' into a datetime format
date_format= "%Y-%m-%dT%H:%M:%S-%f"
date2 = date.apply(lambda x:dt.datetime.strptime(x, date_format))
#print(date2)

#Figure 1, plot time series 
fig1=plt.figure()
plt.plot(date2,pH) #to plot date correctly, must be in datetime formate
plt.xlabel('Date')
plt.ylabel('pH')
plt.title('pH Time Series')
plt.show()

#remove mean from the data
pHmean=np.mean(pH)
#print(pHmean)
nn=len(pH)
pH_nomean = np.zeros(nn)

for i in range (0,nn):
	pH_nomean[i] = pH[i] - pHmean #this single array does not have an index
	#pH_nomean2 = pd.Series(pH_nomean) #pd.Series means it has an index

#Figure 2, plot time series, pH mean removed
fig2=plt.figure()
plt.plot(date2,pH_nomean)
plt.xlabel('Date')
plt.ylabel('pH, mean removed')
plt.title('pH Time Series, Mean Removed')
plt.show()

#Figure 3, plot time series, show trend
fig3=plt.figure()
#fit linear model
X = [i for i in range(0,nn)]
X = np.reshape(X, (len(X), 1))
# print(X)
y = pH_nomean
model = LinearRegression()
model.fit(X, y)
#calculate trend
trend = model.predict(X)
# plot trend
plt.plot(y)
plt.plot(trend)
plt.xlabel('Date')
plt.ylabel('pH, mean removed')
plt.title('pH Time Series, Mean Removed, Trend')
plt.show()

#Figure 4, plot time series, remove trend and the mean
fig4=plt.figure()
detrended = [y[i]-trend[i] for i in range(0,nn)]
plt.plot(detrended)
plt.xlabel('Date')
plt.ylabel('pH, mean and trend removed')
plt.title('pH Time Series, Mean and Trend Removed')
plt.show()

#determine the frequencies for the spectrum
delt=2
T_length=nn
pi=np.pi
omega0=2.*pi/(T_length*delt)
mm=int(nn/2+1)
freq=np.zeros(mm)
#
for i in range(0,mm):
	freq[i] = i*omega0

zz=np.fft.rfft(detrended,n=nn)
fourier_amp=np.sqrt((np.real(zz)**2+np.imag(zz)**2))
#fourier_phase=180.*arctan2(np.imag(zz),np.real(zz))/pi

#plot Fourier amplitudes as a function of frequency
fig5=plt.figure()
plt.loglog(freq,fourier_amp)
plt.xlabel('$\omega$ (radians/???)')
plt.ylabel('Fourier Amplitude')
plt.title('Fourier Amplitude vs. Time')
plt.show()

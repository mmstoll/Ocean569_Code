# plots the digital fourier transform results
# for a 5-day (120 hr) and 90 day (2160 hr)
# digital running mean filter 
# and also a 5-day triangular filter
#
# this will be applied to the MBARI M1 mooring dataset, with
# a length of 156884 points; the filters are padded with zeros to this length
#
import numpy as np
import matplotlib.pyplot as plt
#
path_out='/Users/riser/Desktop/ocean.569A/MBARI.running.mean.and.triangular.filters.jpg'
#
# set the constants
#
nn=156884
mm=int(nn/2+1)
T_length=nn
pi=np.pi
omega0=2.*pi/T_length
#
tt=np.zeros(nn)
yy=np.zeros(nn)
zz=np.zeros(nn)
freq0=np.zeros(mm)
freq1=np.zeros(mm)
#
T_length=nn
pi=np.pi
omega0=2.*pi/T_length
#
# set the filter widths and period marks for the plots
#
n_len=[121,2159,121]
n_plots=len(n_len)
omega_chk=2.*pi*nn*np.array([1./n_len[0],1./n_len[1],1./98.5])
#
# begin the fft and plotting loop
#
fig=plt.figure(figsize=(12,10))
#
for ii in range(0,n_plots):
    plt.subplot(n_plots,1,ii+1)
# 
# scale the plots so that the sum of the weights equals 1
#
    amp=np.array([1./(n_len[ii]+1),1./(n_len[ii]+1),1])    
#
# define the filters
#
    sum_signal=0.
#
    for i in range(0,nn):
      tt[i]=i
      if i <= n_len[ii] and ii < 2:
        signal=amp[ii]
      else:
        signal=0.
      if i <= n_len[ii]//2 and ii == 2:
        signal=i*4./(n_len[ii])**2
      if i > n_len[ii]//2 and i <= n_len[ii] and ii == 2:
        signal=(2./n_len[ii])*(1.-(i-n_len[ii]/2)*4/((n_len[ii])**2))
      if i > n_len[ii] and ii ==2:
        signal=0.
      yy[i]=signal
      sum_signal+=yy[i]
#
# define the frequencies of the transform (+ frequencies only)
#
    for i in range(0,mm):
      freq0[i]=i*omega0
      freq1[i]=freq0[i]*T_length
#
# carry out the FFT
# 
    zz=np.fft.rfft(yy,n=nn)
    zz_real=zz.real
    zz_imag=zz.imag
    zz_mag=np.sqrt(zz_real**2+zz_imag**2)
    if ii == 2:
        zz_mag=zz_mag/sum_signal
    mmaxx=np.max(zz_mag)
#
# plot the 3 filters on semilog plots
#
    fs=12
    plt.xlim([5.,4.e4])
    plt.xlim([5.,1.e5])
    plt.ylim([0,1.2])
    plt.semilogx(freq1,zz_mag,color='purple')
#
# annotate the plots
#
    if ii == 0:
        plt.text(105.,1.08,'120 hr (5 day) running mean filter',fontsize=fs)
        plt.text(0.9e4,0.6,'5-day period',fontsize=0.8*fs)
    if ii == 1:
        plt.text(105.,1.08,'2160 hr (90 day) running mean filter',fontsize=fs)
        plt.text(550.,0.6,'90-day period',fontsize=0.8*fs)
    if ii == 2:
        plt.text(105.,1.08,'120 hr (5 day) triangular filter',fontsize=fs)
        plt.text(1.1e4,0.6,'4.1-day period',fontsize=0.8*fs)
        plt.xlabel('$\omega\it{T}$',fontsize=fs)
        plt.text(8.5e2,-0.27,'[continues to 4.9x10$^5$]',fontsize=fs)
    plt.yticks([0,0.5,1,1.2],['0',' ','1',' '])
    plt.text(2.25,1.1,'|$\it{H}(\omega$)|',fontsize=fs)
    plt.text(10.5,0.4,'$\it{T}$ = 156884 hr = 17.9 yr',fontsize=fs)
    plt.plot([omega_chk[ii],omega_chk[ii]],[0.25,0.88],'--g') 
    plt.grid()
# 
# save the plot if desired
#
#plt.savefig(path_out)
#
plt.show()
#
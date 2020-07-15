"""

Data: Temeprature and Salinity time series from SIO Scripps Pier
	Salinity: measured in PSU at the surface (~0.5m) and at depth (~5m)
	Temp: measured in degrees C at the surface (~0.5m) and at depth (~5m)
- Timestamp included beginning in 1990

"""

# imports 
import numpy as np
import scipy.stats as ss

def var_fft(data_set):
	ll = len(data_set)
	if ll < 2000:
		delt = 30
	if ll >= 2000:
		delt = 1
	data_fft = data_set

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

	return(freq,spec,spec_amp,fft,delt,freq_T,freq_nyquist)

# =================================================================================

def conf_int(spec_data, conf_lim, df):
	spec_data_len = len(spec_data)
	conf_l = np.zeros(spec_data_len)
	conf_h = np.zeros(spec_data_len)
	conf_above = (1.-conf_lim)/2
	conf_below = 1.-(1.-conf_lim)/2
	for i in range (0,len(spec_data)):
		conf_l[i] = spec_data[i] * ss.chi2.ppf([conf_above], df)
		conf_h[i] = spec_data[i] * ss.chi2.ppf([conf_below], df)
	return(conf_l, conf_h)

# =================================================================================

def band_average(fft_var1,fft_var2,frequency,n_av,delt):
	# fft_var1 and fft_var2 are the inputs computed via fft
	# they can be the same variable or different variables
	# n_av is the number of bands to be used for smoothing (nice if it is an odd number)
	# this function is limnited to 100,000 points but can easily be modified
    nmax=100000
    T_length = (len(fft_var1) * 2 - 2)

	# define some variables and arrays
    n_spec=len(fft_var1)
    n_av2=int(n_av//2+1) #number of band averages/2 + 1
    spec_amp_av=np.zeros(nmax)
    spec_phase_av=np.zeros(nmax)
    freq_av=np.zeros(nmax)
	# average the lowest frequency bands first (with half as many points in the average)
    sum_low_amp=0.
    sum_low_phase=0.
    count=0
    spectrum_amp=np.absolute(fft_var1*np.conj(fft_var2))/(2.*np.pi*T_length*delt)
    spectrum_phase=np.angle(fft_var1*np.conj(fft_var2), deg=True) #/(2.*np.pi*T_length*delt),deg=True) #don't know if I need the 2pi/Tdeltt here...
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
        # omega0 = 2.*np.pi/(T_length*delt)
	# contract the arrays
    spec_amp_av=spec_amp_av[0:count]
    spec_phase_av=spec_phase_av[0:count]
    freq_av=freq_av[0:count]
    return spec_amp_av,spec_phase_av,freq_av,count

def band_average_uncorr(fft_var1,fft_var2,frequency,n_av,delt):
	# fft_var1 and fft_var2 are the inputs computed via fft
	# they can be the same variable or different variables
	# n_av is the number of bands to be used for smoothing (nice if it is an odd number)
	# this function is limnited to 100,000 points but can easily be modified
    nmax=100000
    # T_length = (len(fft_var1) * 2 - 2)

	# define some variables and arrays
    n_spec=len(fft_var1)
    n_av2=int(n_av//2+1) #number of band averages/2 + 1
    spec_amp_av=np.zeros(nmax)
    spec_phase_av=np.zeros(nmax)
    freq_av=np.zeros(nmax)
	# average the lowest frequency bands first (with half as many points in the average)
    sum_low_amp=0.
    sum_low_phase=0.
    count=0
    spectrum_amp=np.absolute(fft_var1*np.conj(fft_var2))#/(2.*np.pi*T_length*delt)
    spectrum_phase=np.angle(fft_var1*np.conj(fft_var2), deg=True) #/(2.*np.pi*T_length*delt),deg=True) #don't know if I need the 2pi/Tdeltt here...
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
        # omega0 = 2.*np.pi/(T_length*delt)
	# contract the arrays
    spec_amp_av=spec_amp_av[0:count]
    spec_phase_av=spec_phase_av[0:count]
    freq_av=freq_av[0:count]
    return spec_amp_av,spec_phase_av,freq_av,count



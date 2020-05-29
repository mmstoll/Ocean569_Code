#
# estimate and plot a chi-squared distribution with error bounds
#
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
#
# file_out='/Users/riser/Desktop/ocean.569A/chi-sq.nu=10.confidence.limits.jpg'
#
# define the number of degrees of freedom
#
nu=10
# 
# define the arrays
#
nn=1000
xx=np.zeros(nn)
yy=np.zeros(nn)
xx=np.linspace(-1,30,nn)
yy=stats.chi2.pdf(xx,nu,loc=0, scale=1)
#
# define the confidence limits
#
plt.plot(xx,yy)
plt.show()
# 
conf_lim=0.95
conf_above=(1.-conf_lim)/2
conf_below=1.-(1.-conf_lim)/2.
mark_rt = stats.chi2.ppf(conf_below,nu)
mark_lf = stats.chi2.ppf(conf_above,nu)
print (mark_lf,mark_rt)
#
# set up the plot object
#
fig=plt.figure(figsize=(12,5))
plt.xlim([-0.75,30])
plt.ylim([0,0.16])
plt.plot(xx,yy,'k')
plt.plot([nu,nu],[0,stats.chi2.pdf(nu,nu,loc=0,scale=1)],'--r')
plt.xlabel('$\chi^2$', fontsize=17)
plt.ylabel(r'Probability', fontsize=17)
plt.title(r'$\chi^2\ \mathrm{distribution}, \nu$ = %d' % nu, fontsize=13)
plt.fill_between(xx,0,yy,where=(np.array(xx)>min(xx))&(np.array(xx)<=mark_lf),facecolor='magenta')
plt.fill_between(xx,0,yy,where=(np.array(xx)>mark_lf)&(np.array(xx)< mark_rt),facecolor='lemonchiffon')
plt.fill_between(xx,0,yy,where=(np.array(xx)>mark_rt)&(np.array(xx)<= max(xx)),facecolor='magenta')
#
# annotate the plot
#
fs=10
plt.text(3.5,0.004,'2.5% outlier',fontsize=fs,color='magenta')
plt.text(17.5,0.004,'2.5% outlier',fontsize=fs,color='magenta')
plt.text(0.75,0.04,'$\chi^2_{0.975} = %.2f$'%mark_lf,fontsize=fs)
plt.text(20.7,0.015,'$\chi^2_{0.025} = %.2f$'%mark_rt,fontsize=fs)
plt.text(15.5,0.125,'$\chi^2_{0.975} \leq \chi^2 \leq \chi^2_{0.025}$',fontsize=fs)
plt.text(15.5,0.105,'$3.25 \leq \chi^2 \leq 20.48$',fontsize=fs)
plt.grid(True)
# 
# save the figure and close
#
#plt.savefig(file_out)
plt.show()

# estimate chi-squared upper and lower error bars
# for a given set of confidence limits
#
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
#
# path_out='/Users/riser/Desktop/ocean.569A/chi-sq.log-log.error.bars.jpg'
#
# set up the plot parameters
#
nn=500
xx=np.zeros(nn)
yy_lower0=np.zeros(nn)
yy_upper0=np.zeros(nn)
frac=0.1
#
fig=plt.figure(figsize=(8,5))
plt.xlim([2,50])
plt.ylim([-2,3])
fs=15
#
# define the set of confidence limits and their colors
#
conf_set=[0.95,0.80,0.70]
color_set=['red','darkblue','goldenrod']
#
# loop through the set of confidence limits
# fine the upper and lower bound for each
#
for j in range(0,3):
  for i in range(0,nn-2):
    ii=frac*(i+2)
    nu=ii
    d_mean=nu
    conf_lim=conf_set[j]
    xx[i]=ii
    conf_above=(1.-conf_lim)/2
    conf_below=1.-(1.-conf_lim)/2.
    mark_rt = stats.chi2.ppf(conf_below,nu)
    mark_lf = stats.chi2.ppf(conf_above,nu)
    yy_upper_nolog=mark_rt
    yy_lower_nolog=mark_lf
    yy_upper0[i]=(yy_upper_nolog-d_mean)/d_mean
    yy_lower0[i]=(d_mean-yy_lower_nolog)/d_mean
  ind_set=nn-2
  xx=xx[0:ind_set]
  yy_lower1=yy_lower0[0:ind_set]
  yy_upper1=yy_upper0[0:ind_set]
#
# plot each curve
#
  plt.plot(xx,yy_upper1,color=color_set[j])
  plt.plot(xx,-yy_lower1,color=color_set[j],linestyle=(0,(1,1)))
#
# finish the plot with appropriate labels
#
plt.plot([2,50],[0,0],'k')
plt.text(32,0.7,'Upper error bound',fontsize=fs)
plt.text(32,-0.9,'Lower error bound',fontsize=fs)
plt.xlabel(r'$\nu$',fontsize=fs)
plt.ylabel('Fraction of the mean value',fontsize=fs)
plt.grid()
#
fs=13
plt.plot([21,23],[2.8,2.8],color=color_set[0])
plt.text(23.8,2.75,r'$\alpha$ = 0.95',fontsize=fs)
plt.plot([21,23],[2.5,2.5],color=color_set[1])
plt.text(23.8,2.45,r'$\alpha$ = 0.80',fontsize=fs)
plt.plot([21,23],[2.2,2.2],color=color_set[2])
plt.text(23.8,2.15,r'$\alpha$ = 0.70',fontsize=fs)
#
#plt.savefig(path_out)
plt.show()
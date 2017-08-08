import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,logspace,cov,array,argsort,log,histogram,log10,median,arange,genfromtxt,zeros,percentile
import cPickle as pickle

path_to_files='/Users/dmelgar/USGSFF/param/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
run_regression=True
Nbins=500

#Get the catalog


#Run regression
Niter=100e3
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)
Asr,ksr,mcmcsr=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=Niter,burn=50e3,fix_exponent=False,dependent_variable='slip_rate')
Apr,kpr,mcmcpr=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=Niter,burn=50e3,fix_exponent=False,dependent_variable='potency_rate')










#Make some plots

#Sommerville scaling
Mo=logspace(19,23)
tr=4e-7*(Mo**(1./3))


fig=plt.figure(figsize=(17,7))












ax=fig.add_subplot(121)

trace = [mcmcsr.trace('A')[:],
        mcmcsr.trace('k')[:],
        mcmcsr.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#228B22',marker='o',s=40)




#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#228B22',marker='o',s=40)



plt.ylim([0.001,0.2])
plt.xlim([1e19,1e23])
plt.ylabel('Mean slip rate (m/s)',fontsize=14)

yticks = [0.01,0.1,1]
labels = [0.01,0.1,1]
plt.yticks(yticks, labels)

#add magnitude ticks
ax2 = ax.twiny()
ax2.set_xscale('log')
new_tick_locations = array([10**(1.5*7+9.1),
                            10**(1.5*7.5+9.1),
                            10**(1.5*8.0+9.1),
                            10**(1.5*8.5+9.1),
                            10**(1.5*9.0+9.1)])
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels([7.0,7.5,8.0,8.5,9.0])
ax2.set_xlabel("Moment Magnitude",fontsize=14)
ax2.minorticks_off()







ax=fig.add_subplot(122)

trace = [mcmcpr.trace('A')[:],
        mcmcpr.trace('k')[:],
        mcmcpr.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)


i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#228B22',marker='o',s=40)




#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,tr,'#808080',lw=2)
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_potency_rate[i],c='#228B22',marker='o',s=40)



plt.ylim([1e6,1e10])
plt.xlim([1e19,1e23])
plt.ylabel('Mean max pot. rate (m**3/s)',fontsize=14)

yticks = [1e6,1e7, 1e8, 1e9, 1e10]
labels = yticks
plt.yticks(yticks, labels)

#add magnitude ticks
ax2 = ax.twiny()
ax2.set_xscale('log')
new_tick_locations = array([10**(1.5*7+9.1),
                            10**(1.5*7.5+9.1),
                            10**(1.5*8.0+9.1),
                            10**(1.5*8.5+9.1),
                            10**(1.5*9.0+9.1)])
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(new_tick_locations)
ax2.set_xticklabels([7.0,7.5,8.0,8.5,9.0])
ax2.set_xlabel("Moment Magnitude",fontsize=14)
ax2.minorticks_off()








plt.subplots_adjust(top=0.94,bottom=0.08,hspace=0.25,wspace=0.2)

plt.show()

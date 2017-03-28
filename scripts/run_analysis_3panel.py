import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,logspace,cov,array,argsort,log,histogram,log10,median,arange

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
run_regression=True
Nbins=500

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)

#Run regression
if run_regression:
    Atr,ktr,mcmctr=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='rise_time')
    Apw,kpw,mcmcpw=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='pulse_length')








#Make some plots

#Sommerville scaling
Mo=logspace(19,23)
tr=4e-7*(Mo**(1./3))


fig=plt.figure(figsize=(17,4.5))

ax=fig.add_subplot(131)
ax.set_xscale('log')


i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#228B22',marker='o',s=40)


plt.legend(['interplate','upper','lower','n/a'])


#Replot
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#228B22',marker='o',s=40)



plt.ylim([1,4.5])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Mean rupture speed (km/s)',fontsize=14)

yticks = [1, 2, 3, 4,5]
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










ax=fig.add_subplot(132)

trace = [mcmctr.trace('A')[:],
        mcmctr.trace('k')[:],
        mcmctr.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,tr,'#808080',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#228B22',marker='o',s=40)


plt.legend(['This study','Sommerville 99'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,tr,'#808080',lw=2)
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#228B22',marker='o',s=40)



plt.ylim([0.9,30])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Mean rise time (s)',fontsize=14)

yticks = [1, 2, 5, 10, 20]
labels = [1, 2, 5, 10, 20]
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





ax=fig.add_subplot(133)

trace = [mcmcpw.trace('A')[:],
        mcmcpw.trace('k')[:],
        mcmcpw.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#228B22',marker='o',s=40)


plt.legend(['This study'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#228B22',marker='o',s=40)



plt.ylim([3,80])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Pulse length (km)',fontsize=14)

yticks = [5,10, 20, 50]
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

plt.subplots_adjust(top=0.85,bottom=0.16)

plt.show()

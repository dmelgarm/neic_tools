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
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)

#Run regression
Niter=100e3
Atr,ktr,mcmctr=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=Niter,burn=50e3,fix_exponent=False,dependent_variable='rise_time')
Apw,kpw,mcmcpw=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='pulse_length')
Apr,kpr,mcmcpr=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='slip_rate')
Adu,kdu,mcmcdu=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='duration')
Act,kct,mcmcct=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=50e3,fix_exponent=False,dependent_variable='centroid_time')





#### Get statistics on rise time

path_to_files='/Users/dmelgar/USGSFF/param/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
run_regression=True
Nbins=500

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)

ev=genfromtxt('/Users/dmelgar/USGSFF/Tr_minmax.txt',usecols=0,dtype='S')
tr=genfromtxt('/Users/dmelgar/USGSFF/Tr_minmax.txt')

Tr_min=tr[:,2]+tr[:,3]
Tr_max=tr[:,2]*tr[:,4]+tr[:,3]*tr[:,4] 


#Find rise times, and durations in the order in which they are in Gavin's tr file
tr_min=zeros(neic.Nevents)
tr_max=zeros(neic.Nevents)
duration=zeros(neic.Nevents)
mean_tr=zeros(neic.Nevents)
percent=zeros(neic.Nevents)
rise_times=[]
count=0
for k in range(len(ev)):
    e=ev[k]
    i=where(neic.IDs==e)[0]
    if len(i)>0: #event exists
        i=i[0]
        tr_min[count]=Tr_min[k]
        tr_max[count]=Tr_max[k]
        duration[count]=neic.neic_durations[i]
        mean_tr[count]=neic.mean_rise_times[i]/duration[i]
        percent[count]=percentile(neic.rise_times[i]/Tr_max[k],90)
        rise_times.append(neic.rise_times[i]/duration[count])
        count+=1
    else:
        print 'Event '+e+' not in current dataset'
        
#normalize durations
tr_min=tr_min/duration  
tr_max=tr_max/duration  


############################################







#Make some plots

#Sommerville scaling
Mo=logspace(19,23)
tr=4e-7*(Mo**(1./3))


fig=plt.figure(figsize=(17,11))

ax=fig.add_subplot(331)
ax.set_xscale('log')


i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['interplate','upper','lower','n/a'])


#Replot
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#228B22',marker='o',s=40,edgecolor='k')



plt.ylim([1,4.5])
plt.xlim([1e19,1e23])
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










ax=fig.add_subplot(332)

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
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['This study','Sommerville 99'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,tr,'#808080',lw=2)
i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rise_times[i],c='#228B22',marker='o',s=40,edgecolor='k')



plt.ylim([0.9,30])
plt.xlim([1e19,1e23])
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





ax=fig.add_subplot(333)

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
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['This study'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#228B22',marker='o',s=40,edgecolor='k')



plt.ylim([3,80])
plt.xlim([1e19,1e23])
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








ax=fig.add_subplot(334)

trace = [mcmcdu.trace('A')[:],
        mcmcdu.trace('k')[:],
        mcmcdu.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['This study'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.neic_durations[i],c='#228B22',marker='o',s=40,edgecolor='k')



plt.ylim([6,180])
plt.xlim([1e19,1e23])
plt.ylabel('Source Duration (s)',fontsize=14)

yticks = [10, 20, 50,100]
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
ax2.minorticks_off()










ax=fig.add_subplot(335)

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
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['This study'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_slip_rates[i],c='#228B22',marker='o',s=40,edgecolor='k')



plt.ylim([1e-3,0.2])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel(r'Mean slip rate (m/s)',fontsize=14)

#yticks = [1e7, 1e9, 1e10]
#labels = [r'$10^7'
#plt.yticks(yticks, labels)

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
ax2.minorticks_off()




#Centroid time
#Lingling ye scaling
Mo=logspace(19,23)
Tc=2.58e-6*(Mo**(1./3))

ax=fig.add_subplot(336)

trace = [mcmcct.trace('A')[:],
        mcmcct.trace('k')[:],
        mcmcct.trace('sigma')[:]]

A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)

ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,Tc,'k--',lw=1.5)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#87CEFA',marker='D',s=40,edgecolor='k')

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#DC143C',marker='s',s=40,edgecolor='k')

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#FFD700',marker='^',s=40,edgecolor='k')

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.legend(['This study','Ye et al. 16'],loc=4,frameon=False)


#Replot
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
plt.plot(Mo,Tc,'k--',lw=1.5)

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#87CEFA',marker='D',s=40,edgecolor='k')
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#DC143C',marker='s',s=40,edgecolor='k')
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#FFD700',marker='^',s=40,edgecolor='k')
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.centroid_times[i],c='#228B22',marker='o',s=40,edgecolor='k')


plt.ylim([3,120])
plt.xlim([1e19,1e23])

yticks = [5, 10, 20,50,100]
plt.yticks(yticks, yticks)

plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel(r'Centroid time (s)',fontsize=14)

#yticks = [1e7, 1e9, 1e10]
#labels = [r'$10^7'
#plt.yticks(yticks, labels)

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
ax2.minorticks_off()






ax=fig.add_subplot(337)
ax.hist(mean_tr,facecolor='#3CB371',bins=15)
ax.set_xlabel('Mean rise time / Event duration',fontsize=14)
ax.set_ylabel('Count',fontsize=14)
ax.set_xlim([0,1])




plt.subplots_adjust(top=0.93,bottom=0.08,hspace=0.35,wspace=0.25,left=0.1,right=0.95)

plt.show()

import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,logspace,cov,array,argsort,log,histogram,log10,median,arange
from matplotlib.ticker import MultipleLocator

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
run_regression=True
Nbins=500

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True)







yl=[0,0.3]


fig=plt.figure(figsize=(11,9))

ax=fig.add_axes([0.1,0.1,0.6,0.8])
axbar=fig.add_axes([0.7,0.1,0.15,0.8])


ax.set_xscale('log')

i=where(neic.event_class=='i')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#228B22',marker='o',s=40)


ax.legend(['interplate','upper','lower','n/a'],loc=0)


#Replot
i=where(neic.event_class=='i')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#87CEFA',marker='D',s=40)
i=where(neic.event_class=='u')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_class=='l')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#FFD700',marker='^',s=40)
i=where(neic.event_class=='n/a')[0]
ax.scatter(neic.moments[i],neic.mean_rise_times[i]/neic.event_durations[i],c='#228B22',marker='o',s=40)



ax.set_ylim(yl)
ax.set_xlim([1e19,0.9e23])
ax.set_xlabel('Moment (Nm)',fontsize=14)
ax.set_ylabel(r'$\tau_r/T_d$',fontsize=16)

majorLocator = MultipleLocator(0.1)
minorLocator = MultipleLocator(0.01)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)


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


## bar
freq,bins=histogram(neic.mean_rise_times/neic.event_durations,bins=15)
axbar.barh(bins[0:-1],freq,height=bins[1]-bins[0],facecolor='#FF8C00')
axbar.set_xlim([0,freq.max()+2])
axbar.set_ylim(yl)
axbar.yaxis.set_ticklabels([])
majorLocator = MultipleLocator(5)
axbar.xaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(0.1)
axbar.yaxis.set_major_locator(majorLocator)
axbar.set_xlabel('Event count',fontsize=14)



fig=plt.figure(figsize=(9,9))
ax=fig.add_subplot(111)
ax.set_xscale('log')


plt.scatter(neic.moments,neic.mean_rise_times/neic.event_durations,c='#1E90FF',marker='D',s=40)
Dx=0.05
Dy=0.001
for k in range(len(neic.moments)):
    dx=neic.moments[k]*Dx
    dy=neic.mean_rupture_velocity[k]*Dy
    plt.annotate(xy=(neic.moments[k]+dx,neic.mean_rise_times[k]/neic.event_durations[k]+dy),s=neic.names[k],fontsize=10)


plt.ylim(yl)
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel(r'$\tau_r/T_d$',fontsize=14)

yticks = [0.1,0.2,0.3]
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

plt.show()

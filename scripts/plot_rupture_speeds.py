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
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True)










fig=plt.figure(figsize=(18,9))
ax=fig.add_subplot(121)
ax.set_xscale('log')

#2-sigma region

i=where(neic.event_type=='S')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#1E90FF',marker='D',s=40)

i=where(neic.event_type=='N')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_type=='T')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40)

plt.legend(['Strike-slip','Normal','Thrust'],loc=2)


#replot things after legend
i=where(neic.event_type=='S')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#1E90FF',marker='D',s=40)
i=where(neic.event_type=='N')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_type=='T')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40)

plt.ylim([1.0,4.5])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Mean rupture speed (km/s)',fontsize=14)

yticks = [1,2,3,4]
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




ax=fig.add_subplot(122)
ax.set_xscale('log')

i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#87CEFA',marker='D',s=40)

i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#FFD700',marker='^',s=40)

i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_rupture_velocity[i],c='#228B22',marker='o',s=40)


plt.legend(['interplate','upper','lower','n/a'],loc=2)


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
plt.ylabel('Mean rupture velocity (km/s)',fontsize=14)

yticks = [1,2,3,4]
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





fig=plt.figure(figsize=(9,9))
ax=fig.add_subplot(111)
ax.set_xscale('log')


plt.scatter(neic.moments,neic.mean_rupture_velocity,c='#1E90FF',marker='D',s=40)
Dx=0.05
Dy=0.015
for k in range(len(neic.moments)):
    dx=neic.moments[k]*Dx
    dy=neic.mean_rupture_velocity[k]*Dy
    plt.annotate(xy=(neic.moments[k]+dx,neic.mean_rupture_velocity[k]+dy),s=neic.names[k],fontsize=10)


plt.ylim([1,4.5])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Mean rupture velocity (km/s)',fontsize=14)

yticks = [1,2,3,4]
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

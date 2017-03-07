import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,logspace

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14


#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2)

#Run regression
#A,k=neic.run_regression(inversion_type='linear',fix_exponent=False)
#Ab,kb,mcmc=neic.run_regression(inversion_type='bayesian',fix_exponent=False)

#Make some plots


fig=plt.figure(figsize=(9,9))
ax=fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

i=where(neic.event_type=='S')[0]
#plt.scatter(neic.mean_slip[i],neic.mean_rise_times[i],c='#1E90FF',marker='D',s=40)
plt.scatter(neic.max_slip[i],neic.mean_rise_times[i],c='#1E90FF',marker='D',s=40)

i=where(neic.event_type=='N')[0]
#plt.scatter(neic.mean_slip[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40)
plt.scatter(neic.max_slip[i],neic.mean_rise_times[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_type=='T')[0]
#plt.scatter(neic.mean_slip[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40)
plt.scatter(neic.max_slip[i],neic.mean_rise_times[i],c='#FFD700',marker='^',s=40)


plt.legend(['Strike-slip','Normal','Thrust'],loc=2)

plt.ylim([1,30])
plt.xlim([0.1,65])
plt.xlabel('Mean slip (m)',fontsize=14)
plt.ylabel('Mean rise time (s)',fontsize=14)

yticks = [1, 2, 5, 10, 20]
labels = [1, 2, 5, 10, 20]
plt.yticks(yticks, labels)

plt.show()
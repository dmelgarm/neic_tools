from numpy import genfromtxt,argsort,where,zeros,arange,percentile
import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl


path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
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
        

fig=plt.figure(figsize=(22,4))
ax=fig.add_subplot(111)
ax.boxplot(rise_times,whis=[5,95],sym='')
plt.scatter(arange(len(duration))+1,tr_max,c='b',lw=0,s=25)
ax.set_xlabel('Event number',fontsize=14)
ax.set_ylabel('Normalized rise time',fontsize=14)
ax.set_ylim([0,0.8])
xticks = range(0,150,10)
labels = xticks
plt.xticks(xticks, labels)
plt.subplots_adjust(top=0.98,right=0.98,left=0.05,bottom=0.13)

fig=plt.figure(figsize=(14,4))
ax=fig.add_subplot(131)
ax.hist(mean_tr/tr_max,facecolor='#3CB371',bins=15)
ax.set_xlabel('Mean rise time / Maximum allowed rise time',fontsize=14)
ax.set_ylabel('Count',fontsize=14)
ax.set_xlim([0,1])

ax=fig.add_subplot(132)
ax.hist(percent,facecolor='#FF9933',bins=15)
ax.set_xlabel('90th pctile rise time / Maximum allowed rise time',fontsize=14)
ax.set_ylabel('Count',fontsize=14)
ax.set_xlim([0,1])

ax=fig.add_subplot(133)
ax.hist(mean_tr,facecolor='#FFD700',bins=15)
ax.set_xlabel('Mean rise time / Event duration',fontsize=14)
ax.set_ylabel('Count',fontsize=14)
ax.set_xlim([0,1])

plt.subplots_adjust(bottom=0.15,top=0.97,left=0.1,right=0.97)
    
plt.show()
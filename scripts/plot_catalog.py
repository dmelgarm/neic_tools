import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
from numpy import zeros,histogram
from obspy.core import UTCDateTime

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2)

#get year
year=zeros(neic.Nevents)
for k in range(neic.Nevents):
    date=UTCDateTime(neic.origin_times[k])
    yr=date.year+date.month/12.+date.day/365.
    year[k]=yr
    

plt.figure(figsize=(7,4.5))

# left,bottom,width,heigth
ax1 = plt.axes([0.1,0.12, 0.7,0.6])
ax2 = plt.axes([0.1,0.72, 0.7,0.25])
ax3 = plt.axes([0.8,0.12, 0.15,0.6])

ax=ax1
markerline, stemlines, baseline=ax.stem(year,neic.magnitudes,facecolor='k')
plt.setp(markerline, 'markerfacecolor', '#FF4500','lw',1)
plt.setp(stemlines, 'color','k', 'linewidth', 1)
ax.set_ylim([6.6,9.2])
ax.set_xlim([1989.5,2017.5])
ax.set_xlabel('Year',fontsize=14)
ax.set_ylabel('Magnitude',fontsize=14)

ax=ax2
hist,bins=histogram(year,bins=20)
ax.bar(bins[0:-1],hist,width=bins[1]-bins[0],facecolor='#6495ED')
ax.set_xlim([1989.5,2017.5])
ax.set_ylabel('Event count',fontsize=14)
ax.tick_params(labelbottom='off')



ax=ax3
hist,bins=histogram(neic.magnitudes,bins=20)
ax.barh(bins[0:-1],hist,height=bins[1]-bins[0],facecolor='#3CB371')
ax.set_ylim([6.6,9.2])
ax.set_xlabel('Event count',fontsize=14)
ax.tick_params(labelleft='off')
ax.xaxis.set_ticklabels(['0','','10','','20',''])

plt.show()
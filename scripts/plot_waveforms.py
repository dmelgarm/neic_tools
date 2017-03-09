from obspy import read
from matplotlib import pyplot as plt
import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
import matplotlib as mpl
from numpy import genfromtxt
from matplotlib.ticker import MultipleLocator

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
#Get the catalog
#neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True)

#Sample waveforms
st1=read(u'/Users/dmelgar/Maule2010/GPS/proc/cons.LXE.sac')
st2=read(u'/Users/dmelgar/Iquique2014/GPS/proc/psga.LXE.sac')
st3=read(u'/Users/dmelgar/NewZealand2016/GPS/proc/kaik.LXE.sac')
st4=read(u'/Users/dmelgar/Nepal2015/GPS/PPP/KKN4.LXN.sac')
st5=read(u'/Users/dmelgar/Nicoya2012/GPS/proc/qsec.LXN.sac')

#moment rates
mr1=genfromtxt(u'/Users/dmelgar/USGSFF/stfs/maule.mr')
mr2=genfromtxt(u'/Users/dmelgar/USGSFF/stfs/iquique.mr')
mr3=genfromtxt(u'/Users/dmelgar/USGSFF/stfs/kaikoura.mr')
mr4=genfromtxt(u'/Users/dmelgar/USGSFF/stfs/nepal.mr')
mr5=genfromtxt(u'/Users/dmelgar/USGSFF/stfs/nicoya.mr')

fig=plt.figure(figsize=(18,2.0))


ax=plt.subplot(151)
xl=[-20,200]
ax.fill_between(mr1[:,0], mr1[:,1]/1.e7, facecolor='#606060',alpha=0.5)
ax.set_xlim(xl)
minorLocator = MultipleLocator(5)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.1e20)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_ticks([1e20,2e20,3e20])
ax.yaxis.set_ticklabels(['10','20','30'])
ax.set_ylabel(r'$\dot M\; (10^{19}Nm/s)$',fontsize=14)
####
ax2 = ax.twinx()
dt=23705
yl=[-6,1]
ax2.plot(st1[0].times()-dt,st1[0].data,'k')
ax.plot([0,0],ax.get_ylim(),'--k')
ax2.plot(ax.get_xlim(),[0,0],'--k')
ax2.set_ylim(yl)
ax2.set_xlim(xl)
minorLocator = MultipleLocator(0.2)
ax2.yaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(2)
ax2.yaxis.set_major_locator(majorLocator)
###



ax=plt.subplot(152)
xl=[-20,150]
ax.fill_between(mr2[:,0], mr2[:,1]/1.e7, facecolor='#606060',alpha=0.5)
ax.set_xlim(xl)
minorLocator = MultipleLocator(5)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.2e19)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_ticks([2e19,4e19,6e19,8e19])
ax.yaxis.set_ticklabels(['2','4','6','8'])
####
ax2 = ax.twinx()
dt=85637
yl=[-1.1,0.1]
ax2.plot(st2[0].times()-dt,st2[0].data,'k')
ax.plot([0,0],ax.get_ylim(),'--k')
ax2.plot(ax.get_xlim(),[0,0],'--k')
ax2.set_ylim(yl)
ax2.set_xlim(xl)
minorLocator = MultipleLocator(0.05)
ax2.yaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(0.5)
ax2.yaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(50)
ax2.xaxis.set_major_locator(majorLocator)
###



ax=plt.subplot(153)
xl=[-10,130]
ax.fill_between(mr3[:,0], mr3[:,1]/1.e7, facecolor='#606060',alpha=0.5)
ax.set_xlim(xl)
minorLocator = MultipleLocator(4)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.2e19)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_ticks([1e19,2e19,3e19,4e19])
ax.yaxis.set_ticklabels(['1','2','3','4'])
ax.set_xlabel(r'$Seconds$',fontsize=14)
####
ax2 = ax.twinx()
dt=65
yl=[-0.2,1.15]
ax2.plot(st3[0].times()-dt,st3[0].data,'k')
ax.plot([0,0],ax.get_ylim(),'--k')
ax2.plot(ax.get_xlim(),[0,0],'--k')
ax2.set_ylim(yl)
ax2.set_xlim(xl)
minorLocator = MultipleLocator(0.05)
ax2.yaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(0.5)
ax2.yaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(50)
ax2.xaxis.set_major_locator(majorLocator)
###


ax=plt.subplot(154)
xl=[-10,100]
ax.fill_between(mr4[:,0], mr4[:,1]/1.e7, facecolor='#606060',alpha=0.5)
ax.set_xlim(xl)
minorLocator = MultipleLocator(4)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.1e19)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_ticks([1e19,2e19])
ax.yaxis.set_ticklabels(['1','2'])
####
ax2 = ax.twinx()
dt=4319
yl=[-2.1,0.2]
ax2.plot(st4[0].times()-dt,st4[0].data,'k')
ax.plot([0,0],ax.get_ylim(),'--k')
ax2.plot(ax.get_xlim(),[0,0],'--k')
ax2.set_ylim(yl)
ax2.set_xlim(xl)
minorLocator = MultipleLocator(0.1)
ax2.yaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(0.5)
ax2.yaxis.set_major_locator(majorLocator)
###


ax=plt.subplot(155)
xl=[-5,70]
ax.fill_between(mr5[:,0], mr5[:,1]/1.e7, facecolor='#606060',alpha=0.5)
ax.set_xlim(xl)
minorLocator = MultipleLocator(2)
ax.xaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.1e19)
ax.yaxis.set_minor_locator(minorLocator)
ax.yaxis.set_ticks([1e19])
ax.yaxis.set_ticklabels(['1'])
####
ax2 = ax.twinx()
dt=52952
yl=[-0.5,0.1]
ax2.plot(st5[0].times()-dt,st5[0].data-0.39,'k')
ax.plot([0,0],ax.get_ylim(),'--k')
ax2.plot(ax.get_xlim(),[0,0],'--k')
ax2.set_ylim(yl)
ax2.set_xlim(xl)
minorLocator = MultipleLocator(0.02)
ax2.yaxis.set_minor_locator(minorLocator)
majorLocator = MultipleLocator(0.2)
ax2.yaxis.set_major_locator(majorLocator)
majorLocator = MultipleLocator(20)
ax2.xaxis.set_major_locator(majorLocator)
ax2.set_ylabel(r'$Displ. (m)$',fontsize=14)
###




plt.subplots_adjust(left=0.05,right=0.95,bottom=0.25,top=0.95,wspace=0.4)


plt.show()
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

#Run regression
if run_regression:
    A,k,mcmc=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=1000e3,burn=100e3,fix_exponent=False,dependent_variable='pulse_length')
    print '\n'
    print A
    print k
else:
    A=9.05037048766e-06
    k=0.279444862572







#Make some plots

# plot MCMC results
fig = plt.figure(figsize=(6.5, 6.5))


# left,bottom,width,heigth
ax1 = plt.axes([0.25,0.15, 0.6,0.6])
ax2 = plt.axes([0.25,0.75, 0.6,0.2])
ax3 = plt.axes([0.05,0.15, 0.2,0.6])
cbax= plt.axes([0.28,0.18,0.01,0.4])



#extrac traces for plotting
trace = [mcmc.trace('A')[:],
        mcmc.trace('k')[:],
        mcmc.trace('sigma')[:]]

#plot scatter and contour sigmas

ax=ax1
im=ax.hist2d(trace[0],trace[1],bins=Nbins,cmap=plt.cm.gist_heat_r)
ax.set_xlim([trace[0].min(),trace[0].max()])
ax.set_ylim([trace[1].min(),trace[1].max()])
ax.scatter(trace[0].mean(),trace[1].mean(),marker='s',c='#3CB371',lw=1,s=50)
ax.set_xlabel(r'$log(A)$',fontsize=16)
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
ax.set_ylabel(r'$k$',fontsize=16)
ax.annotate(r'$t_r=A M_0^k$',xy=(-5.42,0.45),fontsize=17)
ax.annotate('freq.',xy=(-8.19,0.424),fontsize=14)
fig.colorbar(im[3], cax=cbax)

ax=ax2
hist,bins=histogram(trace[0],bins=Nbins)
ax.yaxis.tick_right()
ax.fill_between(bins[:-1],hist,facecolor='#B0C4DE',alpha=0.4)
ax.plot(bins[:-1],hist,lw=1.0,color='k')
ax.set_ylabel('Frequency',fontsize=14)
ax.yaxis.set_label_position("right")
ax.tick_params(labelbottom='off')
ax.set_xlim([trace[0].min(),trace[0].max()])
ax.plot([trace[0].mean(),trace[0].mean()],ax.get_ylim(),lw=2,c='#3CB371')
ax.yaxis.set_ticks(arange(0,6001,1000))
ax.yaxis.set_ticklabels(['','1000','','3000','','5000',''])


ax=ax3
hist,bins=histogram(trace[1],bins=Nbins)
ax.fill_betweenx(bins[:-1],hist,facecolor='#B0C4DE',alpha=0.4)
ax.plot(hist,bins[:-1],color='k')
ax.set_xlabel('Frequency',fontsize=14)
ax.tick_params(labelleft='off')
ax.plot(ax.get_xlim(),[trace[1].mean(),trace[1].mean()],lw=2,c='#3CB371')
ax.set_ylim([trace[1].min(),trace[1].max()])
ax.invert_xaxis()
ax.xaxis.set_ticks(arange(0,6001,1000))
ax.xaxis.set_ticklabels(['','','2000','','4000','','6000'])
for label in ax.get_xmajorticklabels():
    label.set_rotation(50)
    label.set_horizontalalignment("center")







# Plot data and regression
Mo=logspace(19,23)

fig=plt.figure(figsize=(18,9))
ax=fig.add_subplot(121)
ax.set_yscale('log')
ax.set_xscale('log')

#2-sigma region
A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mu = yfit.mean(0)
sig = 2 * yfit.std(0)
plt.plot(Mo,10**mu,'k',lw=2)


i=where(neic.event_type=='S')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#1E90FF',marker='D',s=40)

i=where(neic.event_type=='N')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40)

i=where(neic.event_type=='T')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40)

plt.legend(['This study','Strike-slip','Normal','Thrust'],loc=2)


#replot things after legend
ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)
i=where(neic.event_type=='S')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#1E90FF',marker='D',s=40)
i=where(neic.event_type=='N')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#DC143C',marker='s',s=40)
i=where(neic.event_type=='T')[0]
plt.scatter(neic.moments[i],neic.mean_pulse_lengths[i],c='#FFD700',marker='^',s=40)

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




ax=fig.add_subplot(122)
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


plt.legend(['This study','interplate','upper','lower','n/a'],loc=2)


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





fig=plt.figure(figsize=(9,9))
ax=fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

ax.fill_between(Mo, 10**(mu - sig), 10**(mu + sig), color='lightgray')
plt.plot(Mo,10**mu,'k',lw=2)


plt.scatter(neic.moments,neic.mean_pulse_lengths,c='#1E90FF',marker='D',s=40)
Dx=0.05
Dy=0.05
for k in range(len(neic.moments)):
    dx=neic.moments[k]*Dx
    dy=neic.mean_pulse_lengths[k]*Dy
    plt.annotate(xy=(neic.moments[k]+dx,neic.mean_pulse_lengths[k]+dy),s=neic.names[k],fontsize=10)


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

plt.show()

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
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2)

#Run regression
if run_regression:
    Ai,ki,mcmci=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=100e3,burn=20e3,fix_exponent=False,select_event_type='i',dependent_variable='slip')
    Au,ku,mcmcu=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=100e3,burn=20e3,fix_exponent=False,select_event_type='u',dependent_variable='slip')
    Al,kl,mcmcl=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=100e3,burn=20e3,fix_exponent=False,select_event_type='l',dependent_variable='slip')
    Ana,kna,mcmcna=neic.run_regression(inversion_type='bayesian',prior='uninformative',Niter=100e3,burn=20e3,fix_exponent=False,select_event_type='n/a',dependent_variable='slip')
else:
    A=9.05037048766e-06
    k=0.279444862572
    
    
#Sommerville scaling
Mo=logspace(19,23)
    
    
#get 2 sigma regions
trace = [mcmci.trace('A')[:],
        mcmci.trace('k')[:],
        mcmci.trace('sigma')[:]]
A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mui = yfit.mean(0)
sigi = 2 * yfit.std(0)

trace = [mcmcu.trace('A')[:],
        mcmcu.trace('k')[:],
        mcmcu.trace('sigma')[:]]
A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
muu = yfit.mean(0)
sigu = 2 * yfit.std(0)

trace = [mcmcl.trace('A')[:],
        mcmcl.trace('k')[:],
        mcmcl.trace('sigma')[:]]
A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
mul = yfit.mean(0)
sigl = 2 * yfit.std(0)

trace = [mcmcna.trace('A')[:],
        mcmcna.trace('k')[:],
        mcmcna.trace('sigma')[:]]
A, k = trace[:2]
yfit = A[:, None] + k[:, None] * log10(Mo)
muna = yfit.mean(0)
signa = 2 * yfit.std(0)
  
      
    
fig=plt.figure(figsize=(9,9))
  
ax=fig.add_subplot(111)
ax.set_yscale('log')
ax.set_xscale('log')

plt.plot(Mo,10**mui,'#000080',lw=2)
plt.plot(Mo,10**mul,'#FFD700',lw=2)
plt.plot(Mo,10**muu,'#DC143C',lw=2)
plt.plot(Mo,10**muna,'#228B22',lw=2)
plt.legend(['interplate','lower','upper','n/a'],loc=2)



#Replot
ax.fill_between(Mo, 10**(mui - sigi), 10**(mui + sigi), color='#000080',alpha=0.5)
plt.plot(Mo,10**mui,'#000080',lw=2)
ax.fill_between(Mo, 10**(mul - sigl), 10**(mul + sigl), color='#FFD700',alpha=0.3)
plt.plot(Mo,10**mul,'#FFD700',lw=2)
ax.fill_between(Mo, 10**(muu - sigu), 10**(muu + sigu), color='#DC143C',alpha=0.2)
plt.plot(Mo,10**muu,'#DC143C',lw=2)
ax.fill_between(Mo, 10**(muna - signa), 10**(muna + signa), color='#228B22',alpha=0.3)
plt.plot(Mo,10**muna,'#228B22',lw=2)



i=where(neic.event_class=='i')[0]
plt.scatter(neic.moments[i],neic.mean_slip[i],c='#000080',marker='D',s=40,zorder=90)
i=where(neic.event_class=='l')[0]
plt.scatter(neic.moments[i],neic.mean_slip[i],c='#FFD700',marker='^',s=40,zorder=90)
i=where(neic.event_class=='u')[0]
plt.scatter(neic.moments[i],neic.mean_slip[i],c='#DC143C',marker='s',s=40,zorder=90)
i=where(neic.event_class=='n/a')[0]
plt.scatter(neic.moments[i],neic.mean_slip[i],c='#228B22',marker='o',s=40,zorder=90)



plt.ylim([0.2,40])
plt.xlim([1e19,1e23])
plt.xlabel('Moment (Nm)',fontsize=14)
plt.ylabel('Mean slip (m)',fontsize=14)

yticks = [0.2,0.5,1.0,2.0,5.0,10.0,20.0]
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
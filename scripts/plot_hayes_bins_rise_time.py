import sys
sys.path.append('/Users/dmelgar/code/python/neic_tools')
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,array,r_,zeros,mean,arange,meshgrid
from matplotlib.ticker import MultipleLocator
from mudpy import forward
from mtspec import mtspec,wigner_ville_spectrum
from scipy.integrate import trapz
from matplotlib import colors
import matplotlib.cm as cmx
from matplotlib.ticker import MultipleLocator
from matplotlib import rcParams

path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
run_regression=True
Nbins=500

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True)

i=where(neic.magnitudes<7.5)[0]
rt1=neic.mean_rise_times[i].mean()
ms1=neic.mean_slip[i].mean()
msr1=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.1) & (neic.magnitudes<7.6))[0]
rt2=neic.mean_rise_times[i].mean()
ms2=neic.mean_slip[i].mean()
msr2=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.2) & (neic.magnitudes<7.7))[0]
rt3=neic.mean_rise_times[i].mean()
ms3=neic.mean_slip[i].mean()
msr3=neic.mean_slip_rates[i].mean()


i=where((neic.magnitudes>7.3) & (neic.magnitudes<7.8))[0]
rt4=neic.mean_rise_times[i].mean()
ms4=neic.mean_slip[i].mean()
msr4=neic.mean_slip_rates[i].mean()


i=where((neic.magnitudes>7.4) & (neic.magnitudes<7.9))[0]
rt5=neic.mean_rise_times[i].mean()
ms5=neic.mean_slip[i].mean()
msr5=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.5) & (neic.magnitudes<8.0))[0]
rt6=neic.mean_rise_times[i].mean()
ms6=neic.mean_slip[i].mean()
msr6=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.6) & (neic.magnitudes<8.1))[0]
rt7=neic.mean_rise_times[i].mean()
ms7=neic.mean_slip[i].mean()
msr7=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.7) & (neic.magnitudes<8.2))[0]
rt8=neic.mean_rise_times[i].mean()
ms8=neic.mean_slip[i].mean()
msr8=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.8) & (neic.magnitudes<8.3))[0]
rt9=neic.mean_rise_times[i].mean()
ms9=neic.mean_slip[i].mean()
msr9=neic.mean_slip_rates[i].mean()

i=where((neic.magnitudes>7.9) & (neic.magnitudes<8.4))[0]
rt10=neic.mean_rise_times[i].mean()
ms10=neic.mean_slip[i].mean()
msr10=neic.mean_slip_rates[i].mean()

rt=[rt1,rt2,rt3,rt4,rt5,rt6,rt7,rt8,rt9,rt10]
ms=[ms1,ms2,ms3,ms4,ms5,ms6,ms7,ms8,ms9,ms10]
msr=[msr1,msr2,msr3,msr4,msr5,msr6,msr7,msr8,msr9,msr10]
M=arange(7.25,8.16,0.1)

fig=plt.figure(figsize=(18,6))
ax=fig.add_subplot(131)
plt.errorbar(M,array(rt),xerr=0.25,c='k')
plt.scatter(M,array(rt),c='r',marker='o',s=20)
plt.xlabel('Magnitude bin',fontsize=14)
plt.ylabel('Mean rise time',fontsize=14)

ax=fig.add_subplot(132)
plt.errorbar(M,array(ms),xerr=0.25,c='k')
plt.scatter(M,array(ms),c='r',marker='o',s=20)
plt.xlabel('Magnitude bin',fontsize=14)
plt.ylabel('Mean slip',fontsize=14)

ax=fig.add_subplot(133)
plt.errorbar(M,array(msr),xerr=0.25,c='k')
plt.scatter(M,array(msr),c='r',marker='o',s=20)
plt.xlabel('Magnitude bin',fontsize=14)
plt.ylabel('Mean peak slip rate',fontsize=14)
ax.set_ylim([0,2])





stf_type='dreger'
zeta=0.25
#zeta=0.2
rate=5
t,stf1=forward.build_source_time_function(rt1,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf2=forward.build_source_time_function(rt2,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf3=forward.build_source_time_function(rt3,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf4=forward.build_source_time_function(rt4,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf5=forward.build_source_time_function(rt5,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf6=forward.build_source_time_function(rt6,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf7=forward.build_source_time_function(rt7,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf8=forward.build_source_time_function(rt8,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf9=forward.build_source_time_function(rt9,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)
t,stf10=forward.build_source_time_function(rt10,0.02,30,stf_type=stf_type,dreger_falloff_rate=rate,zeta=zeta)



#rescale to mean slip
stf1=stf1*(ms1/trapz(stf1,t))
stf2=stf2*(ms2/trapz(stf2,t))
stf3=stf3*(ms3/trapz(stf3,t))
stf4=stf4*(ms4/trapz(stf4,t))
stf5=stf5*(ms5/trapz(stf5,t))
stf6=stf6*(ms6/trapz(stf6,t))
stf7=stf7*(ms7/trapz(stf7,t))
stf8=stf8*(ms8/trapz(stf8,t))
stf9=stf9*(ms9/trapz(stf9,t))
stf10=stf10*(ms10/trapz(stf10,t))

#Get mean slip rate
msr_theoretical=[stf1.max(),stf2.max(),stf3.max(),stf4.max(),stf5.max(),stf6.max(),stf7.max(),
                    stf8.max(),stf9.max(),stf10.max()]
                    
#plot
plt.figure()
plt.scatter(msr,msr_theoretical)
plt.legend(['Mean peak slip rate (m/s)'])
plt.plot([1,2],[1,2])
plt.xlabel('From NEIC subfault data',fontsize=14)
plt.ylabel('From theoretical slip rate function',fontsize=14)                
                    

#zero pad and other business
st1=r_[zeros(500),stf1]
st2=r_[zeros(500),stf2]
st3=r_[zeros(500),stf3]
st4=r_[zeros(500),stf4]
st5=r_[zeros(500),stf5]
st6=r_[zeros(500),stf6]
st7=r_[zeros(500),stf7]
st8=r_[zeros(500),stf8]
st9=r_[zeros(500),stf9]
st10=r_[zeros(500),stf10]



#get spec
s1, f, jackknife, _, _ = mtspec(data=st1, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s2, f, jackknife, _, _ = mtspec(data=st2, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s3, f, jackknife, _, _ = mtspec(data=st3, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s4, f, jackknife, _, _ = mtspec(data=st4, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s5, f, jackknife, _, _ = mtspec(data=st5, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s6, f, jackknife, _, _ = mtspec(data=st6, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s7, f, jackknife, _, _ = mtspec(data=st7, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s8, f, jackknife, _, _ = mtspec(data=st8, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s9, f, jackknife, _, _ = mtspec(data=st9, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)
s10, f, jackknife, _, _ = mtspec(data=st10, delta=0.02, time_bandwidth=3.5,number_of_tapers=5, statistics=True)


#Spectrograms
#smoothing_filter='gauss',filter_width=10
wv1 = wigner_ville_spectrum(st1, delta=0.02, time_bandwidth=3.5,smoothing_filter='gauss',filter_width=50)
wv10 = wigner_ville_spectrum(st10, delta=0.02, time_bandwidth=3.5,smoothing_filter='gauss',filter_width=50)
wv1=abs(wv1)
wv10=abs(wv10)
twv=arange(0,len(st1)*0.02-0.01,0.02)-500*0.02
T,F=meshgrid(twv,f[::-1])


# Plot the source time functions
# Using contourf to provide my colorbar info, then clearing the figure
Z = [[0,0],[0,0]]
cm = plt.get_cmap('brg') 
levels = arange(7.25,8.25+0.01,0.01)
plt.figure()
c = plt.contourf(Z, levels, cmap=cm)
plt.clf()
cNorm  = colors.Normalize(vmin=7.25, vmax=8.25)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)



rcParams['xtick.major.size'] = 6.0
rcParams['xtick.major.width'] = 0.5
rcParams['xtick.minor.size'] = 4
rcParams['xtick.minor.width'] = 0.5
rcParams['ytick.major.size'] = 6.0
rcParams['ytick.major.width'] = 0.5
rcParams['ytick.minor.size'] = 4
rcParams['ytick.minor.width'] = 0.5


fig=plt.figure(figsize=(16,6))
ax=fig.add_subplot(121)
colorVal1 = scalarMap.to_rgba(7.25)
colorVal2 = scalarMap.to_rgba(7.35)
colorVal3 = scalarMap.to_rgba(7.45)
colorVal4 = scalarMap.to_rgba(7.55)
colorVal5 = scalarMap.to_rgba(7.65)
colorVal6 = scalarMap.to_rgba(7.75)
colorVal7 = scalarMap.to_rgba(7.85)
colorVal8 = scalarMap.to_rgba(7.95)
colorVal9 = scalarMap.to_rgba(8.05)
colorVal10 = scalarMap.to_rgba(8.15)
ax.plot(t,stf1,color=colorVal1,lw=2)
ax.plot(t,stf2,color=colorVal2,lw=2)
ax.plot(t,stf3,color=colorVal3,lw=2)
ax.plot(t,stf4,color=colorVal4,lw=2)
ax.plot(t,stf5,color=colorVal5,lw=2)
ax.plot(t,stf6,color=colorVal6,lw=2)
ax.plot(t,stf7,color=colorVal7,lw=2)
ax.plot(t,stf8,color=colorVal8,lw=2)
ax.plot(t,stf9,color=colorVal9,lw=2)
ax.plot(t,stf10,color=colorVal10,lw=2)
ax.set_xlim([0,12])
ax.set_xlabel(r'$Time (s)$',fontsize=16)
ax.set_ylabel(r'$Slip\;\;rate (m/s)$',fontsize=16)

xmajorLocator = MultipleLocator(2)
xminorLocator = MultipleLocator(0.2)
ymajorLocator = MultipleLocator(0.2)
yminorLocator = MultipleLocator(0.05)
ax.xaxis.set_major_locator(xmajorLocator)
ax.yaxis.set_major_locator(ymajorLocator)
ax.xaxis.set_minor_locator(xminorLocator)
ax.yaxis.set_minor_locator(yminorLocator)

#cb=ax.colorbar(c)
#cb.set_label('Bin center magnitude')

ax=fig.add_subplot(122)
ax.loglog(f,s1,color=colorVal1)
ax.loglog(f,s2,color=colorVal2)
ax.loglog(f,s3,color=colorVal3)
ax.loglog(f,s4,color=colorVal4)
ax.loglog(f,s5,color=colorVal5)
ax.loglog(f,s6,color=colorVal6)
ax.loglog(f,s7,color=colorVal7)
ax.loglog(f,s8,color=colorVal8)
ax.loglog(f,s9,color=colorVal9)
ax.loglog(f,s10,color=colorVal10)
ax.loglog([1./13,1./13],[1e-9,1e2],'--',c='#606060',lw=2)
ax.set_xlabel(r'$Frequency (Hz)$',fontsize=16)
ax.set_ylabel(r'$Power (m^2/s^2/Hz)$',fontsize=16)
ax.set_xlim([0.03,20])
ax.set_ylim([5e-7,9e-1])
cb=fig.colorbar(c)
cb.set_label(r'$Bin center magnitude$',fontsize=16)


#plot spectrograms
vmin=min(wv1.min(),wv10.min())
vmax=max(wv1.max(),wv10.max())
cmap=plt.cm.jet

fig=plt.figure(figsize=(12,6))
ax=fig.add_subplot(211)
ax.set_yscale('log')
cb=ax.pcolormesh(T,F,wv1,vmin=vmin,vmax=vmax,cmap=cmap)
ax.set_ylim([0.03,20])
ax.set_ylabel('Frequency',fontsize=14)
cbar=fig.colorbar(cb)
cbar.set_label('Wigner-Ville amplitude')
ax.annotate(r'$7.0\leq M_w<7.5$',xy=(20,4.5),color='w',fontsize=18)



ax=fig.add_subplot(212)
ax.set_yscale('log')
cb=ax.pcolormesh(T,F,wv10,vmin=vmin,vmax=vmax,cmap=cmap)
ax.set_ylim([0.03,20])
cbar=fig.colorbar(cb)
ax.set_xlabel('Time (s)',fontsize=14)
ax.set_ylabel('Frequency',fontsize=14)
cbar.set_label('Wigner-Ville amplitude')
ax.annotate(r'$7.9\leq M_w<8.4$',xy=(20,4.5),color='w',fontsize=18)

plt.show()





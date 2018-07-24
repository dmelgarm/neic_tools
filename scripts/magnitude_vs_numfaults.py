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


#Run regression
Niter=100e3
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)


for k in range(len(neic.magnitudes)):
    plt.scatter(neic.magnitudes[k],len(neic.segment_data[k])
    
plt.show()  
from numpy import genfromtxt,argsort,where,zeros,arange,percentile
import sys
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl


path_to_files='/Users/dmelgar/USGSFF/param/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 18
mpl.rcParams['ytick.labelsize'] = 18
run_regression=True
Nbins=500

#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2,get_stfs=True,percent_cutoff_vrup=0.3)


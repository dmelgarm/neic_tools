# -*- coding: utf-8 -*-
from pymc import Model, Normal, HalfNormal
import pymc as pm 
from neic_tools import neic_catalog
from matplotlib import pyplot as plt
import matplotlib as mpl
from numpy import where,logspace,log10,array
import corner


path_to_files='/Users/dmelgar/Downloads/PARAM_FILES/'
catalog_file='/Users/dmelgar/USGSFF/catalog.txt'
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14


#Get the catalog
neic=neic_catalog(path_to_files,catalog_file,percent_cutoff=0.2)

def myModel(x, y_observe): 
    # Priors for unknown model parameters
    #Note tau is the 1/sigma^2, e.g 1/σ^2
    A_prior = pm.distributions.Normal('A', mu = -5, tau = 1)
    B_prior = pm.distributions.Normal('B', mu = 1./3, tau = 1)
    
    #You can try different distributions on the variance of the noise
    #sigma = basic_model.HalfNormal('sigma', 1)
    #tau_prior = pm.distributions.Gamma("tau", alpha=0.001, beta=0.001)
    #Here I use sigma = 0.1
    tau_prior = 1.0/1.0**2
    
    # Expected value of outcome
    @pm.deterministic(plot=False)
    def f(x = x, A=A_prior, B=B_prior):
        return A+B*x
    
    # Likelihood (sampling distribution) of observations
    Y_obs = pm.distributions.Normal('Y_obs', mu = f, tau = tau_prior, value = y_observe, observed=True)
    return locals()

model = pm.Model(myModel(log10(neic.moments),log10(neic.mean_rise_times)))

#Here I also calculate MAP to get the best solution
model = pm.Model(myModel(log10(neic.moments),log10(neic.mean_rise_times)))
map_ = pm.MAP(model)
map_.fit()
#print the MAP solution
print map_.A_prior.value, map_.B_prior.value

#Use MCMC to get the distribution of parameters
mcmc = pm.MCMC(myModel(log10(neic.moments),log10(neic.mean_rise_times)))
mcmc.sample(iter=100000, burn=20000)
pm.Matplot.plot(mcmc)

samples = array([mcmc.A_prior.trace(),mcmc.B_prior.trace()]).T
tmp = corner.corner(samples[:,:], labels=['A','B'],bins=50)

plt.show()
#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, print_function
import emcee
import numpy as np
import matplotlib.pyplot as pl
from scipy import stats

#np.random.seed(123)

# Define the model.
def model(p):
    m, b = p
    return lambda x0: m * x0 + b

def log_prior(theta):
    m, b, s = theta
    if s < 0:
        return -np.inf  # log(0)
    else:
        #return -1.5 * np.log(1 + m ** 2) - np.log(s)
        return 0

def log_likelihood(theta, x, y, xerr, yerr):
    m, b, s = theta
    model = m*x + b
    sigma2 = s**2 + yerr**2 + m**2*xerr**2

    return -0.5 * np.sum(np.log(2 *np.pi*sigma2) + (y-model)**2 / sigma2)

def log_probfn(theta, x, y, xerr, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return log_prior(theta) + log_likelihood(theta, x, y, xerr, yerr)

# Generate data.
truth = [1.5, 4.]
N = 30
x_true = 50 * np.random.rand(N)
y_true = model(truth)(x_true)

# "Observe" the data.
xerr, yerr = 2.0, 10.0

x_obs = stats.norm(x_true, xerr).rvs(N)
y_obs = stats.norm(y_true, yerr).rvs(N)

# Set up the sampler.
nwalkers, ndim = 100, 3
p0 = np.random.random((nwalkers, ndim))
#p0 = np.append(truth, x_obs)
#p0 = [p0 + np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probfn,
                                args=[x_obs, y_obs, xerr, yerr])

# Burn in.
print("Burning in.")
p0, lnprob0, state = sampler.run_mcmc(p0, 500)
sampler.reset()

# Sample.
print("Sampling.")
sampler.run_mcmc(p0, 1000)

# Print results.
samples = sampler.flatchain
print("m = {0} ± {1}".format(np.mean(samples[:, 0]), np.std(samples[:, 0])))
print("b = {0} ± {1}".format(np.mean(samples[:, 1]), np.std(samples[:, 1])))
print("s = {0} ± {1}".format(np.mean(samples[:, 2]), np.std(samples[:, 2])))

# Plot results.
fig = pl.figure(figsize=(10, 5))

xl = np.linspace(0,50)

# plot crediable region
m, b = samples.T[:2]
yfit = m[:,None]* xl + b[:,None] # This creates individual points
mu = yfit.mean(0)
sig = 2 * yfit.std(0) # 2sigma confidence
pl.fill_between(xl, mu - sig, mu + sig, color='lightgray')
pl.plot(xl, mu, '-k', label='Fit')

# plot true line
pl.plot(xl, truth[0]*xl+truth[1], color="r", lw=2, alpha=0.8, label='True')

# and points
pl.errorbar(x_true, y_obs, xerr=xerr, yerr=yerr, fmt='o', label='Obs')
pl.show()

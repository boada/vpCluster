from __future__ import division, print_function
import pylab as pyl
import h5py as hdf
import emcee

# Define the model.
def model(p):
    m, b = p
    return lambda x0: m * x0 + b

def log_prior(theta):
    m, b, s = theta
    if s < 0:
        return -pyl.inf  # log(0)
    else:
        #return -1.5 * np.log(1 + m ** 2) - np.log(s)
        return 0

def log_likelihood(theta, x, y, xerr, yerr):
    m, b, s = theta
    model = m*x + b
    sigma2 = s**2 + yerr**2 + m**2*xerr**2

    return -0.5 * pyl.sum(pyl.log(2 *pyl.pi*sigma2) + (y-model)**2 / sigma2)

def log_probfn(theta, x, y, xerr, yerr):
    lp = log_prior(theta)
    if not pyl.isfinite(lp):
        return -pyl.inf
    return log_prior(theta) + log_likelihood(theta, x, y, xerr, yerr)

def mkError(x, x_err):
    return 0.434 * x_err/x

with hdf.File('./../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value
    data.sort(order='ID')

richnessData = pyl.genfromtxt('./boada_rich.txt', names=True, dtype=None)
richnessData.sort(order='name')

x_obs = pyl.log10(richnessData['lambda'])
xerr = mkError(richnessData['lambda'], richnessData['lambda_err'])
y_obs = data['MASS']
yerr = data['MASS_err']

f = pyl.figure(figsize=(7,7*(pyl.sqrt(5.)-1.0)/2.0))
ax = f.add_subplot(111)

# low mass
eb = pyl.errorbar(x_obs[:2], y_obs[:2], xerr=xerr[:2], yerr=yerr[:2],
        fmt='o',mfc='lightgray', mec='#e24a33', color='#e24a33', markersize=10)
eb[-1][0].set_linestyle('--')
eb[-1][1].set_linestyle('--')

# high mass
pyl.errorbar(x_obs[2:], y_obs[2:], xerr=xerr[2:], yerr=yerr[2:],
        fmt='o', color='#e24a33', markersize=10, label='This Work')

ax.set_xlabel('Log Richness')
ax.set_ylabel(r'Log Cluster Mass $(M_{\odot})$')

# add the Rozo2010 points
rozo = pyl.genfromtxt('./rozo2010_points', names=True, dtype=None)
ax.errorbar(pyl.log10(rozo['lambda']), pyl.log10(rozo['mass']*1e15),
        xerr=mkError(rozo['lambda'], rozo['lambda_err']), fmt='o',
        label='Rozo+2010')

# add the masses from Sifon2015
sifonMasses = pyl.array([14.6e14, 4.5e14,17.0e14, 7.0e14])
sifonMass_err = pyl.array([3.1e14, 1.9e14, 2.8e14, 2.0e14])
sifonRichness = pyl.array([190.996, 135.59, 202.636, 185.764])
ax.errorbar(pyl.log10(sifonRichness), pyl.log10(sifonMasses),
        yerr=mkError(sifonMasses, sifonMass_err), fmt='o', color='#7a68a6',
        label='Sifon+2015')

# add the ML predicted points
with hdf.File('./../MLmethods/ML_predicted_masses_shifty.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    ML = dset.value
    ML.sort(order='ID')

y_obs = ML['ML_pred']
yerr = ML['ML_pred_err']
#ax.errorbar(x_obs[2:], y_obs_ML[2:], xerr=xerr[2:], yerr=yerr_ML[2:],
#        fmt='o', color='#467821', markersize=10, label='ML$_{\sigma, N_{gal}}$')

#eb = pyl.errorbar(x_obs[:2], y_obs_ML[:2], xerr=xerr[:2], yerr=yerr_ML[:2],
#        fmt='o',mfc='lightgray', mec='#467821', color='#467821', markersize=10)
eb[-1][0].set_linestyle('--')
eb[-1][1].set_linestyle('--')


# Set up the sampler.
nwalkers, ndim = 100, 3
p0 = pyl.random((nwalkers, ndim))
sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probfn,
        args=[x_obs[2:], y_obs[2:], xerr[2:], yerr[2:]])

# Burn in.
print("Burning in.")
p0, lnprob0, state = sampler.run_mcmc(p0, 500)
sampler.reset()

# Sample.
print("Sampling.")
sampler.run_mcmc(p0, 1000)

# Print results.
samples = sampler.flatchain
print("m = {0} +/- {1}".format(pyl.mean(samples[:, 0]), pyl.std(samples[:, 0])))
print("b = {0} +/- {1}".format(pyl.mean(samples[:, 1]), pyl.std(samples[:, 1])))
print("s = {0} +/- {1}".format(pyl.mean(samples[:, 2]), pyl.std(samples[:, 2])))

xl = pyl.linspace(0.8,2.4)
# plot crediable region
m, b = samples.T[:2]
yfit = m[:,None]* xl + b[:,None] # This creates individual points
mu = yfit.mean(0)
sig = yfit.std(0) # 2sigma confidence
pyl.fill_between(xl, mu - sig, mu + sig, color='lightgray', zorder=0)
pyl.plot(xl, mu, '-k', label='Fit', zorder=0)

# add the Rykoff2012 relation
x = pyl.logspace(0.8,2.4)
lny = 1.48 + 1.06*pyl.log(x/60.)
y = pyl.exp(lny) * 1e14/0.7
ax.plot(pyl.log10(x),pyl.log10(y), ls=':', c='k', zorder=0, label='Rykoff+2012')

# add the Farahi2016 relation
lny = 0.44 + 1.31* pyl.log(x/30.)
y = pyl.exp(lny) * 1e14
ax.plot(pyl.log10(x),pyl.log10(y), ls='--', c='k', zorder=0, label='Farahi+2016')


pyl.legend(loc='lower right', ncol=2)
pyl.show()


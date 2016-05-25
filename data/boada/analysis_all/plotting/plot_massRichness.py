from __future__ import division
import pylab as pyl
import h5py as hdf

def mkError(x, x_err):
    return 0.434 * x_err/x

with hdf.File('./../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value
    data.sort(order='ID')

richnessData = pyl.genfromtxt('./boada_rich.txt', names=True, dtype=None)
richnessData.sort(order='name')


mass = data['MASS']
yerr = data['MASS_err']

f = pyl.figure(figsize=(5,5*(pyl.sqrt(5.)-1.0)/2.0))
ax = f.add_subplot(111)

ax.errorbar(pyl.log10(richnessData['lambda']), mass,
        xerr=mkError(richnessData['lambda'], richnessData['lambda_err']),
        yerr=yerr, fmt='o', color='#e24a33', markersize=10, label='This Work')

ax.set_xlabel('Log Richness')
ax.set_ylabel(r'Log Cluster Mass $[M_{\odot}]$')

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

pyl.legend(loc='lower right')

pyl.show()


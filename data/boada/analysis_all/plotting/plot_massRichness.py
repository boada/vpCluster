from __future__ import division
import pylab as pyl
import h5py as hdf

def mkError(x, x_err):
    return 0.434 * x_err/x

with hdf.File('./../MLmethods/ML_predicted_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    data = dset.value
    data.sort(order='ID')

richnessData = pyl.genfromtxt('./boada_rich.txt', names=True, dtype=None)
richnessData.sort(order='name')


mass = data['ML_pred_3d']
yerr = data['ML_pred_3d_err']
#yerr = mkError(data['MASS'], data['MASS_err'])

ax = pyl.subplot()

ax.errorbar(pyl.log10(richnessData['lambda']), mass,
        xerr=mkError(richnessData['lambda'], richnessData['lambda_err']),
        yerr=yerr, fmt='o', label='This Work')

ax.set_xlabel('Log Richness')
ax.set_ylabel(r'Log Cluster Mass $[M_{\odot}]$')

# add the Rozo2010 points
#rozo = pyl.genfromtxt('./rozo2010_points', names=True, dtype=None)
#ax.errorbar(pyl.log10(rozo['lambda']), pyl.log10(rozo['mass']*1e15),
#        xerr=mkError(rozo['lambda'], rozo['lambda_err']), fmt='o',
#        label='Rozo+2010')

# add the masses from Sifon2015
#sifonMasses = pyl.array([14.6e14, 4.5e14,17.0e14, 7.0e14])
#sifonMass_err = pyl.array([3.1e14, 1.9e14, 2.8e14, 2.0e14])
#sifonRichness = pyl.array([190.996, 135.59, 202.636, 185.764])
#ax.errorbar(pyl.log10(sifonRichness), pyl.log10(sifonMasses),
#        yerr=mkError(sifonMasses, sifonMass_err), fmt='o', color='0.8',
#        label='Sifon+2015')

pyl.legend(loc='lower right')

pyl.show()


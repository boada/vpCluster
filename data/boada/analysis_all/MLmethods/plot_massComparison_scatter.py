from glob import glob
import pylab as pyl
import h5py as hdf

files = glob('ML_predicted_masses*')

# get the power law masses
with hdf.File('../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    results = dset.value

for f,c in zip(files, 'rgb'):

    with hdf.File(f, 'r') as f1:
        dset = f1[f1.keys()[0]]
        ML = dset.value

    pyl.errorbar(results['MASS'], ML['ML_pred'], xerr=results['MASS_err'],
            yerr=ML['ML_pred_err'], fmt='o', color=c, markersize=10,
            label=f.rstrip('.hdf5'))

pyl.xlabel('M$_{pred, PL}$')
pyl.ylabel('M$_{pred, ML}$')
pyl.plot([12,16],[12,16], c='k', zorder=0)



from glob import glob
import pylab as pyl
import h5py as hdf

files = glob('ML_predicted_masses*')

# get the power law masses
with hdf.File('../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    results = dset.value

# make a figure
f = pyl.figure(figsize=(6, 6 * (pyl.sqrt(5.) - 1.0) / 2.0))
ax = f.add_subplot(111)

i = 0
for f, c, l in zip(files, ['#7a68a6', '#e24a33'],
                   ['$ML_{\sigma, N_{gals}}$, Flat HMF',
                    '$ML_{\sigma, N_{gals}}$']):

    if i == 0:
        i += 1
        continue

    with hdf.File(f, 'r') as f1:
        dset = f1[f1.keys()[0]]
        ML = dset.value

    ax.errorbar(results['MASS'],
                ML['ML_pred'],
                xerr=results['MASS_err'],
                yerr=ML['ML_pred_err'],
                fmt='o',
                color=c,
                markersize=10,
                label=l)  #f.rstrip('.hdf5'))

ax.set_xlabel('Log M$_{pred, PL}$')
ax.set_ylabel('Log M$_{pred, ML}$')
ax.plot([12.5, 15.5], [12.5, 15.5], c='k', zorder=0)

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper left')

pyl.show()

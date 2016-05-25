from __future__ import print_function
import h5py as hdf

with hdf.File('./../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    results = dset.value

with hdf.File('./../MLmethods/ML_predicted_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    ML = dset.value


for i in range(results.size):
    print('%s & %d & %d (%d) & %d & %.4f$\pm{%.3f}$ & %d$\pm{%d}$ & '\
    r'%.2f$\pm{%.2f}$ & %.2f$\pm{%.2f}$ \\' % (results['ID'][i],
        results['SOURCES'][i],
    results['Q0'][i],
    results['Q1'][i],
    results['MEMBERS'][i],
    results['Zc'][i],
    results['Zc_err'][i],
    results['LOSVD'][i],
    results['LOSVD_err'][i],
    results['MASS'][i],
    results['MASS_err'][i],
    ML['ML_pred_3d'][i],
    ML['ML_pred_3d_err'][i]
    ))







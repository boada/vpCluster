from __future__ import print_function
import numpy as np
import h5py as hdf

with hdf.File('./../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    results = dset.value

for i in range(results.size):
    print('%s & %d & %d (%d) & %d & %.4f\err{%.3f}{%.3f} & %d\err{%d}{%d} & '\
    r'%.2f\err{%.2f}{%.2f} \\' % (results['ID'][i], results['SOURCES'][i],
    results['Q0'][i], results['Q1'][i], results['MEMBERS'][i],
    results['Zc'][i], results['Zc_err'][i][0], results['Zc_err'][i][1],
    results['LOSVD'][i], results['LOSVD_err'][i][0],
    results['LOSVD_err'][i][1], results['MASS'][i]/1e14,
    results['MASS_err'][i][0]/1e14, results['MASS_err'][i][1]/1e14))







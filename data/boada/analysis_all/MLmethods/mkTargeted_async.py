from multiprocessing import Pool
import h5py as hdf
import numpy as np
from data_handler import mkTruth, mkHalo, mkQs
from halo_handler import find_indices, find_indices_multi
from calc_cluster_props import (updateArray, findClusterRedshift, findLOSV,
                                findLOSVDmcmc, calc_mass_Evrard)
import os

class AsyncFactory:
    def __init__(self, func, cb_func):
        self.func = func
        self.cb_func = cb_func
        self.pool = Pool(maxtasksperchild=10)

    def call(self,*args, **kwargs):
        self.pool.apply_async(self.func, args, kwargs, self.cb_func)

    def wait(self):
        self.pool.close()
        self.pool.join()

def worker(pos, data, center):
    #print "PID: %d \t Value: %d" % (os.getpid(), pos)
    data = updateArray(data)
    data = findClusterRedshift(data)
    #data = findSeperationSpatial(data, center)
    data = findLOSV(data)
    data, sigma_dist = findLOSVDmcmc(data)
    data = calc_mass_Evrard(data, A1D = 1177, alpha = 0.364)
    return pos, data, sigma_dist

def cb_func((pos, data, sigma_dist)):
    if pos % 1000 == 0:
        print pos
    results['IDX'][pos] = pos
    results['CLUSZ'][pos] = data['CLUSZ'][0]
    results['LOSVD'][pos] = data['LOSVD'][0]
    results['MASS'][pos] = data['MASS'][0]
    results['LOSVD_err'][pos] = data['LOSVD_err'][0]
    results['LOSVD_dist'][pos] = sigma_dist[:,0]

if __name__ == "__main__":

    async_worker = AsyncFactory(worker, cb_func)
    halo = mkHalo()
    truth = mkTruth()
    Qs = mkQs()

    # this is the part that makes it realistic or not
    mask = (Qs == 0) | (Qs == 1)
    truth = truth[mask]

    # there are no clusters with mass < 2e11 and more than 5 galaxies
    mask = (halo['upid'] == -1) & (halo['m200c'] > 2e11)
    maskedHalo = halo[mask]
    hids, uniqueIdx = np.unique(maskedHalo['id'], return_index=True)

    # now we find all member halos of the most massive halos
    subHalos = find_indices(halo['upid'], hids)

    # now find all the corresponding galaxies
    gals = find_indices_multi(truth['HALOID'], halo['id'], subHalos)

    # make the results container
    x = [i for i,g in enumerate(gals) if g.size >=5]
    # make the results container
    results = np.zeros((len(x),), dtype=[('IDX', '>i4'),
        ('HALOID', '>i8'),
        ('ZSPEC', '>f4'),
        ('M200c', '>f4'),
        ('CLUSZ', '>f4'),
        ('LOSVD', '>f4'),
        ('MASS', '>f4'),
        ('NGAL', '>i4'),
        ('LOSVD_err', '>f4', (2,)),
        ('LOSVD_dist', '>f4', (10000,))])

    print('do work', len(x), 'clusters to go!')
    keepBad = False
    for j,i in enumerate(x):
        center = (maskedHalo['ra'][uniqueIdx[i]],
                maskedHalo['dec'][uniqueIdx[i]])
        if gals[i].size >= 5:
            async_worker.call(j, truth[gals[i]], center)
            # update results array
            results['HALOID'][j] = maskedHalo['id'][uniqueIdx[i]]
            results['NGAL'][j] = gals[i].size
            results['ZSPEC'][j] = maskedHalo['zspec'][uniqueIdx[i]]
            results['M200c'][j] = maskedHalo['m200c'][uniqueIdx[i]]/0.70
        elif keepBad:
            results['HALOID'][j] = maskedHalo['id'][uniqueIdx[i]]
            results['NGAL'][j] = gals[i].size
            results['ZSPEC'][j] = maskedHalo['zspec'][uniqueIdx[i]]
            results['M200c'][j] = maskedHalo['m200c'][uniqueIdx[i]]/0.70

    async_worker.wait()

    try:
        #os.remove('result_targetedPerfect.hdf5')
        os.remove('result_targetedRealistic.hdf5')
    except OSError:
        pass
    #with hdf.File('result_targetedPerfect.hdf5', 'w') as f:
    #    f['result_targetedPerfect'] = results
    with hdf.File('result_targetedRealistic.hdf5', 'w') as f:
        f['result_targetedRealistic'] = results

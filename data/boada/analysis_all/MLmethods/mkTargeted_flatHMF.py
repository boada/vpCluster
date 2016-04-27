from multiprocessing import Pool
import h5py as hdf
import numpy as np
from data_handler import mkTruth, mkHalo, mkQs
from halo_handler import find_indices, find_indices_multi
from calc_cluster_props import (updateArray, findClusterRedshift, findLOSV,
                                findLOSVDmcmc, calc_mass_Evrard)
import os
import random

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

def randomR():
    ''' randomly choose a location on the unit hemisphere, r = [x,y,z] with
    -1<x<1, -1<y<1, and 0<z<1 and x^2+y^2+z^2 = 1.

    '''

    phi = 2*np.pi*random.random()
    theta = np.arccos(np.random.random())
    r = [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]
    return r

def worker(pos, data, center, rot=False):
    #print "PID: %d \t Value: %d" % (os.getpid(), pos)
    data = updateArray(data)
    data = findClusterRedshift(data)
    #data = findSeperationSpatial(data, center)

    # do the rotations or not
    if rot:
        r = randomR()
        v = np.column_stack([data['VX'], data['VY'], data['VZ']])
        rot = [np.dot(r,vi) for vi in v]
        data['LOSV'] = rot
    else:
        data = findLOSV(data)

    data, sigma_dist = findLOSVDmcmc(data)
    data = calc_mass_Evrard(data, A1D = 1177, alpha = 0.364)
    return pos, data, sigma_dist

def cb_func((pos, data, sigma_dist)):
    if pos % 200 == 0:
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

    # figure out which clusters have an appropriate size
    x = np.array([i for i,g in enumerate(gals) if g.size >=5])

    #### Here is where we are going to make a flat HMF ####
    massBins = np.arange(11.5, 16, 0.5)
    masses = [np.log10(maskedHalo['m200c'][uniqueIdx[i]]/0.70) for i in x]
    indexes = np.digitize(masses, massBins)
    # now we choose some of the indexes uniformly -- How many do we want?
    indxChoice = np.random.randint(1, massBins.size, 1000)

    # figure out where they fit, and pick one, also uniformly
    picks = [np.random.choice(np.where(indexes==i)[0]) for i in indxChoice]

    # make the results container
    results = np.zeros((x[picks].size,), dtype=[('IDX', '>i4'),
        ('HALOID', '>i8'),
        ('ZSPEC', '>f4'),
        ('M200c', '>f4'),
        ('CLUSZ', '>f4'),
        ('LOSVD', '>f4'),
        ('MASS', '>f4'),
        ('NGAL', '>i4'),
        ('LOSVD_err', '>f4', (2,)),
        ('LOSVD_dist', '>f4', (10000,))])

    print('do work', x[picks].size, 'clusters to go!')
    for j,i in enumerate(x[picks]):
        center = (maskedHalo['ra'][uniqueIdx[i]],
                maskedHalo['dec'][uniqueIdx[i]])

        # start the rotations
        #if j > len(x):
        #    async_worker.call(j, truth[gals[i]], center, rot=True)
        #else:
        async_worker.call(j, truth[gals[i]], center, rot=True)
        #worker(j, truth[gals[i]], center, rot=False)

        # update results array
        results['HALOID'][j] = maskedHalo['id'][uniqueIdx[i]]
        results['NGAL'][j] = gals[i].size
        results['ZSPEC'][j] = maskedHalo['zspec'][uniqueIdx[i]]
        # millenium cosmology
        #results['M200c'][j] = maskedHalo['m200c'][uniqueIdx[i]]/0.73
        results['M200c'][j] = maskedHalo['m200c'][uniqueIdx[i]]/0.70

    async_worker.wait()

    print('results')
    try:
        os.remove('surveyComplete_noRotations.hdf5')
    except OSError:
        pass
    with hdf.File('surveyComplete_noRotations.hdf5', 'w') as f:
        f['surveyComplete'] = results

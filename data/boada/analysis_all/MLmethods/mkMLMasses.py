from __future__ import print_function
import pandas as pd
from glob import glob
import numpy as np
from itertools import combinations
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing
import h5py as hdf

### Training Data ###
#####################
with hdf.File('./buzzard_targetedRealistic.hdf5', 'r') as f:
    dset  = f[f.keys()[0]]
    #data = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
    #    'LOSVD_err', 'MASS', 'LOSVD_dist']
    data = dset['ZSPEC', 'M200c', 'LOSVD', 'NGAL']

# You have to clean the data here. This is almost certainly from the fact
# that some of the HALOIDS are repeated at different redshifts. I have a
# prior on the LOSVD calculation which will limit the LOSVD to a maxium.
# Because the clusters are so far apart the LOSVD is super high.

mask = ((np.log10(data['LOSVD']) > 3.12 ) & (data['M200c'] < 10**14.5) |
    (data['LOSVD'] < 50))
maskedDataT = data[~mask]
badData = data[mask]


X = np.column_stack([data['ZSPEC'], data['NGAL'], np.log10(data['LOSVD'])])
y = np.log10(data['M200c'])


#X = preprocessing.scale(X)
X_train, X_test, y_train, y_test = train_test_split(X,y,
        test_size=0.30)

rf = RandomForestRegressor(n_estimators=1000, min_samples_leaf=1,
                                    verbose=1, n_jobs=-1)

rf.fit(X_train, y_train)
mrf = rf.predict(X_test)




for j in range(5):
    print(j)
    magDict = {}
    with hdf.File('./truth/truth'+str(j).zfill(2)+'_Oii.hdf5', 'r') as f:
        dset = f['truth%s_Oii' % (str(j).zfill(2))]
        magDict['u'] = dset['OMAG'][:,0] # u band
        magDict['g'] = dset['OMAG'][:,1] # g band
        magDict['r'] = dset['OMAG'][:,2] # r band
        magDict['i'] = dset['OMAG'][:,3] # i band
        magDict['z'] = dset['OMAG'][:,4] # z band

    # we only want the g mag < 22 galaxies
    mask = np.where(magDict['g'] < 22)[0]
    print(mask.size)
    # create a data array for everything to fit into
    #features = -np.ones((mask.size, 25))
    data = -np.ones((mask.size, 25))

    # mags
    for i, m in enumerate('ugriz'):
        data[:,i] = magDict[m][mask]
    
    # colors
    colors = combinations('ugriz', 2)
    for i, c in enumerate(colors):
        data[:,i+5] = magDict[c[0]][mask] - magDict[c[1]][mask]
        # colors squared
        data[:,i+15] = (magDict[c[0]][mask] - magDict[c[1]][mask])**2

    
    # data = scaler.transform(data)
    data = preprocessing.scale(data)
    data = rfecv.transform(data)
    # now we make the predictions based on the new features we've created
    Qs = clf.predict(data)

    print(np.where(Qs == 0)[0].size/float(Qs.size))
    print(np.where(Qs == 1)[0].size/float(Qs.size))
    print(np.where(Qs == 2)[0].size/float(Qs.size))
        
    # with hdf.File('./truth/truth'+str(j).zfill(2)+'_Oii.hdf5', 'a') as f:
    #     values = -np.ones(magDict['u'].size)
    #     values[mask] = Qs
    #     f['Q'] = values

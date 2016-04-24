from __future__ import print_function
import numpy as np
from itertools import combinations
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn import preprocessing
import h5py as hdf

def pred_ints(model, X, mrf, percentile=68):
    ''' Calculates the prediction intervals of the estimators. '''

    err = []
    for x in range(len(X)):
        preds = []
        for pred in model.estimators_:
            try:
                preds.append(pred.predict(X[x][:,np.newaxis]))
            except ValueError:
                preds.append(pred.predict(X[x].reshape(1,-1)))
        err.append(np.std(preds))
#        err_down.append(np.percentile(preds, (100 - percentile) / 2. ))
#        err_up.append(np.percentile(preds, 100 - (100 - percentile) / 2.))

    err_down = mrf - err
    err_up = mrf + err

    return err_down, err_up

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
#X_train, X_test, y_train, y_test = train_test_split(X,y,
#        test_size=0.30)

rf = RandomForestRegressor(n_estimators=1000, min_samples_leaf=1,
                                    verbose=1, n_jobs=-1)

rf.fit(X, y)
#mrf = rf.predict(X_test)


### our data ###
################
with hdf.File('./../results_cluster.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    obs = dset.value


X_pred = np.column_stack((obs['Zc'], obs['MEMBERS'],
    np.log10(obs['LOSVD'])))

mrf = rf.predict(X_pred)






import numpy as np
import h5py as hdf
from sklearn.ensemble import RandomForestRegressor
#from sklearn.cross_validation import train_test_split
from numpy.lib import recfunctions as rfns
from itertools import permutations, izip
import multiprocessing

def child_initializer(_rf):
    print 'Starting', multiprocessing.current_process().name
    global model
    model = _rf

def updateArray(data):
    ''' Adds the results containers to the data product. '''

    newData = np.zeros(data.size)
    data = rfns.append_fields(data, ['ML_pred_1d', 'ML_pred_2d',
        'ML_pred_2d2', 'ML_pred_3d', 'ML_pred_1d_err', 'ML_pred_2d_err',
        'ML_pred_2d2_err', 'ML_pred_3d_err' ], [newData, newData, newData,
            newData, newData, newData, newData, newData], dtypes='>f4',
            usemask=False)

    return data

def addMasses(results, train, test):
    ''' This does all of the heavy lifting to get the new masses assigned to
    the right places.

    '''

    rf = RandomForestRegressor(n_estimators=1000, min_samples_leaf=1,
            verbose=1, n_jobs=4)
    X = np.log10(train['M200c'])

############
#### 1d ####
############
    y = np.column_stack([np.log10(train['LOSVD'])])
    rf.fit(y, X)
    obs = np.column_stack([np.log10(test['LOSVD'])])
    mrf = rf.predict(obs)

    results['ML_pred_1d'][test['IDX']] = mrf

    # errors
    print('Calculating Error')
    # p = multiprocessing.Pool(maxtasksperchild=1000,
    #         initializer=child_initializer, initargs=([rf]))
    # result = p.map(mp_worker_wrapper, izip(obs, mrf))
    # p.close()
    # p.join()
    # results['ML_pred_1d_err'][test['IDX']] = result

#############
#### 2d #####
#############
    y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC']])
    rf.fit(y, X)
    obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC']])
    mrf = rf.predict(obs)

    results['ML_pred_2d'][test['IDX']] = mrf
    # errors
    print('Calculating Error, 2d')
    # p = multiprocessing.Pool(maxtasksperchild=1000,
    #         initializer=child_initializer, initargs=([rf]))
    # result = p.map(mp_worker_wrapper, izip(obs, mrf))
    # p.close()
    # p.join()
    # results['ML_pred_2d_err'][test['IDX']] = result

#############
#### 2d2 ####
#############
    try:
        y = np.column_stack([np.log10(train['LOSVD']), train['NMEM']])
    except ValueError:
        y = np.column_stack([np.log10(train['LOSVD']), train['NGAL']])
    rf.fit(y, X)
    try:
        obs = np.column_stack([np.log10(test['LOSVD']), test['NMEM']])
    except ValueError:
        obs = np.column_stack([np.log10(test['LOSVD']), test['NGAL']])
    mrf = rf.predict(obs)

    results['ML_pred_2d2'][test['IDX']] = mrf
    # errors
    print('Calculating Error, 2d2')
    # p = multiprocessing.Pool(maxtasksperchild=1000,
    #         initializer=child_initializer, initargs=([rf]))
    # result = p.map(mp_worker_wrapper, izip(obs, mrf))
    # p.close()
    # p.join()
    # results['ML_pred_2d2_err'][test['IDX']] = result

##############
##### 3d #####
##############
    try:
        y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC'],
            train['NMEM']])
    except ValueError:
        y = np.column_stack([np.log10(train['LOSVD']), train['ZSPEC'],
            train['NGAL']])
    rf.fit(y, X)
    try:
        obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC'],
            test['NMEM']])
    except ValueError:
        obs = np.column_stack([np.log10(test['LOSVD']), test['ZSPEC'],
            test['NGAL']])
    mrf = rf.predict(obs)

    results['ML_pred_3d'][test['IDX']] = mrf
    # errors
    print('Calculating Error, 3d')
    # p = multiprocessing.Pool(maxtasksperchild=1000,
    #         initializer=child_initializer, initargs=([rf]))
    # result = p.map(mp_worker_wrapper, izip(obs, mrf))
    # p.close()
    # p.join()
    # results['ML_pred_3d_err'][test['IDX']] = result

    return results

#def mp_pred_ints(model, obs, mrf):
def mp_pred_ints(obs, mrf):
    preds = []
    for pred in model.estimators_:
        try:
            preds.append(pred.predict(obs[:,np.newaxis]))
        except ValueError:
            preds.append(pred.predict(obs.reshape(1,-1)))

    #err_down = mrf - np.std(preds)
    #err_up = mrf + np.std(preds)

    # Bessel corrected std
    err = np.std(preds, ddof=1)

    return err


def mp_worker_wrapper(args):
    return mp_pred_ints(*args)

if __name__ == "__main__":

    ### Train ###
    #############
    with hdf.File('./buzzard_targetedRealistic_flatHMF_shifty.hdf5', 'r') as f:
        dset  = f[f.keys()[0]]
        try:
            train = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
            'LOSVD_err', 'MASS', 'NMEM']
        except ValueError:
            train = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
            'LOSVD_err', 'MASS']

    ### Test ###
    #############
    with hdf.File('./buzzard_targetedRealistic_shifty.hdf5', 'r') as f:
        dset  = f[f.keys()[0]]
        try:
            test = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
            'LOSVD_err', 'MASS', 'NMEM']
        except ValueError:
            test = dset['IDX', 'HALOID', 'ZSPEC', 'M200c', 'NGAL', 'LOSVD',
            'LOSVD_err', 'MASS']

    # add the extra fields
    result = updateArray(test)

    # You have to clean the data here. This is almost certainly from the fact
    # that some of the HALOIDS are repeated at different redshifts. I have a
    # prior on the LOSVD calculation which will limit the LOSVD to a maxium.
    # Because the clusters are so far apart the LOSVD is super high.
    try:
        mask = ((np.log10(train['LOSVD']) > 3.12 ) &
                (train['M200c'] < 10**14.5) | (train['LOSVD'] < 50) |
                (train['NMEM'] < 5))
                print 'Using NMEM'
    except ValueError:
        mask = ((np.log10(train['LOSVD']) > 3.12 ) &
                (train['M200c'] < 10**14.5) | (train['LOSVD'] < 50))
    train = train[~mask]

    try:
        mask = ((np.log10(test['LOSVD']) > 3.12 ) &
                (test['M200c'] < 10**14.5) | (test['LOSVD'] < 50) |
                (test['NMEM'] < 5))
    except ValueError:
        mask = ((np.log10(test['LOSVD']) > 3.12 ) &
                (test['M200c'] < 10**14.5) | (test['LOSVD'] < 50))
    test = test[~mask]

    data = addMasses(result, train, test)
    with hdf.File('buzzard_targetedRealistic_comparison.hdf5', 'w') as f:
        f['predicted masses'] = data
        f.flush()

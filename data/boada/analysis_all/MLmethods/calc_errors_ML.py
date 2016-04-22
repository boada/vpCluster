import numpy as np
import h5py as hdf
from scipy import stats

def error(true, pred, mu):
    ''' Unused, but kept to see how I did it when I wasn't using the Scipy
    functions. Calculates the error on the mean.

    '''
    print true.size,
    if true.size > 1:
        var = np.sum((pred - true - mu)**2) /(true.size - 1)
        sem = np.sqrt(var/true.size)
        return sem
    elif true.size == 1:
        return 0
    else:
        return np.nan

def bias(true, pred):
    ''' unused, but calculates the mean bias. '''

    if true.size > 0:
       return np.sum(pred - true) /true.size
        #return np.median(true)
    else:
        return np.nan

def runningStatistic(stat, true, pred, **kwargs):
    ''' b = bias and s = uncertainty on that bias '''

    bins = np.arange(11.5,16,0.5)
    indx = np.digitize(true, bins)-1
    binNumber = len(bins)

    runningb = []
    runnings = []
    for k in xrange(binNumber):
        print true[indx==k].size,
        b = np.mean(pred[indx==k] - true[indx==k])
        s = stats.sem(pred[indx==k] - true[indx==k])
        print '$%.2f\pm{%.2f}$ &' % (b,s)
        try:
            mean, var, std = stats.mvsdist(pred[indx==k] - true[indx==k])
            #print '$%.2f\pm{%.2f}$ &' % (std.mean(),std.std()),
        except ValueError:
            pass
            #print '$%.2f\pm{%.2f}$ &' % (np.nan,np.nan),
        runningb.append(b)
        runnings.append(s)
    print ''
    return

### Targeted ###
################
with hdf.File('./buzzard_targetedRealistic_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']
# filter bad values
mask = (target['ML_pred_1d'] != 0)
target = target[mask]

for d in [target]:

### Full survey ###
    mean, var, std = stats.mvsdist(np.log10(d['MASS']) - np.log10(d['M200c']))
    s = stats.sem(np.log10(d['MASS']) - np.log10(d['M200c']))
    #print '$%.2f\pm{%.3f}$' % (mean.mean(),s)
    print '$%.2f\pm{%.3f}$' % (std.mean(), std.std())


    print('power law')
    running = runningStatistic(bias, np.log10(d['M200c']),
            np.log10(d['MASS']))
############
#### 1d ####
############
    print('1d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_1d'])

#############
#### 2d #####
#############
    print('2d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_2d'])

##############
##### 3d #####
##############
    print('3d')
    running = runningStatistic(bias, np.log10(d['M200c']),
            d['ML_pred_3d'])

    print '-----'

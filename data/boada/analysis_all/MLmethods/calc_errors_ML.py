import numpy as np
import h5py as hdf
from scipy import stats


def error(true, pred, mu):
    ''' Unused, but kept to see how I did it when I wasn't using the Scipy
    functions. Calculates the error on the mean.

    '''
    print true.size,
    if true.size > 1:
        var = np.sum((pred - true - mu)**2) / (true.size - 1)
        sem = np.sqrt(var / true.size)
        return sem
    elif true.size == 1:
        return 0
    else:
        return np.nan


def bias(true, pred):
    ''' unused, but calculates the mean bias. '''

    if true.size > 0:
        return np.sum(pred - true) / true.size
        #return np.median(true)
    else:
        return np.nan


def runningStatistic(true, pred, **kwargs):
    ''' b = bias and s = scatter '''

    bins = np.arange(11.5, 16, 0.1)
    indx = np.digitize(true, bins)

    runningb = []
    runnings = []
    for k in xrange(1, bins.size):
        print true[indx == k].size, bins[k],
        try:
            mean, var, std = stats.mvsdist(pred[indx == k] - true[indx == k])
            #print '$%.2f\pm{%.2f}$ &' % (std.mean(),std.std()),
            print '$%.3f\pm{%.3f}$ &' % (mean.mean(), mean.std()),
            runnings.append(std.mean())
            runningb.append(mean.mean())
        except ValueError:
            runnings.append(0.0)
            runningb.append(0.0)
            print '$%.3f\pm{%.3f}$ &' % (0.0, 0.0),
        print ''
    return runningb, runnings

### Targeted ###
################
with hdf.File('./buzzard_targetedRealistic_shifty_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_2d2',
                  'ML_pred_3d', 'NMEM']
# filter bad values
mask = (target['ML_pred_1d'] != 0) & (target['NMEM'] >= 5)
target = target[mask]

# make a container for the biases 
# bins = np.arange(11.5,16,0.1)
# results = np.zeros((bins.size,), dtype= [('bins', '>f4'),
#     ('powerlaw_bias', '>f4'),
#     ('ML_bias_1d', '>f4'),
#     ('ML_bias_2d', '>f4'),
#     ('ML_bias_2d2', '>f4'),
#     ('ML_bias_3d', '>f4'),
#     ('powerlaw_scatter', '>f4'),
#     ('ML_scatter_1d', '>f4'),
#     ('ML_scatter_2d', '>f4'),
#     ('ML_scatter_2d2', '>f4'),
#     ('ML_scatter_3d', '>f4')])

# make a container for the biases
bins = np.arange(11.5, 16, 0.1)
results = np.zeros(
    (bins.size - 1, ),
    dtype=[('centers', '>f4'), ('powerlaw_bias', '>f4'), ('ML_bias_1d', '>f4'),
           ('ML_bias_2d', '>f4'), ('ML_bias_2d2', '>f4'),
           ('ML_bias_3d', '>f4'), ('powerlaw_scatter', '>f4'),
           ('ML_scatter_1d', '>f4'), ('ML_scatter_2d', '>f4'),
           ('ML_scatter_2d2', '>f4'), ('ML_scatter_3d', '>f4')])

results['centers'] = (bins[:-1] + bins[1:]) / 2.

for d in [target]:

    ### Full survey ###
    mean, var, std = stats.mvsdist(np.log10(d['MASS']) - np.log10(d['M200c']))
    s = stats.sem(np.log10(d['MASS']) - np.log10(d['M200c']))
    #print '$%.2f\pm{%.3f}$' % (mean.mean(),s)
    print '$%.2f\pm{%.3f}$' % (std.mean(), std.std())

    print('power law')
    running = runningStatistic(np.log10(d['M200c']), np.log10(d['MASS']))
    results['powerlaw_bias'] = running[0]
    results['powerlaw_scatter'] = running[1]

    ############
    #### 1d ####
    ############
    print('1d')
    running = runningStatistic(np.log10(d['M200c']), d['ML_pred_1d'])
    results['ML_bias_1d'] = running[0]
    results['ML_scatter_1d'] = running[1]

    #############
    #### 2d #####
    #############
    print('2d')
    running = runningStatistic(np.log10(d['M200c']), d['ML_pred_2d'])
    results['ML_bias_2d'] = running[0]
    results['ML_scatter_2d'] = running[1]

    #############
    #### 2d2 #####
    #############
    print('2d2')
    running = runningStatistic(np.log10(d['M200c']), d['ML_pred_2d2'])
    results['ML_bias_2d2'] = running[0]
    results['ML_scatter_2d2'] = running[1]

    ##############
    ##### 3d #####
    ##############
    print('3d')
    running = runningStatistic(np.log10(d['M200c']), d['ML_pred_3d'])
    results['ML_bias_3d'] = running[0]
    results['ML_scatter_3d'] = running[1]

    print '-----'

#with hdf.File('biasandScatter.hdf5', 'w') as f:
#    f['biasandscatter'] = results

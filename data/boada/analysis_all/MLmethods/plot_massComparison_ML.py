import pylab as pyl
from astLib import astStats
from sklearn.metrics import median_absolute_error, mean_squared_error
import h5py as hdf
from matplotlib.ticker import AutoMinorLocator


def calc_err(pred, true):
    return (pred - true) / true


golden_mean = (pyl.sqrt(5.) - 1.0) / 2.0
f = pyl.figure(figsize=(10, 10 * golden_mean))

ax1 = pyl.subplot2grid((3, 4), (0, 0), rowspan=2)
ax2 = pyl.subplot2grid((3, 4), (0, 1), rowspan=2, sharex=ax1)
ax3 = pyl.subplot2grid((3, 4), (0, 2), rowspan=2, sharex=ax1, sharey=ax2)
ax4 = pyl.subplot2grid((3, 4), (0, 3), rowspan=2, sharex=ax1, sharey=ax2)

# now for the bottom bits
ax1s = pyl.subplot2grid((3, 4), (2, 0))
ax2s = pyl.subplot2grid((3, 4), (2, 1), sharex=ax1s)
ax3s = pyl.subplot2grid((3, 4), (2, 2), sharex=ax1s, sharey=ax2s)
ax4s = pyl.subplot2grid((3, 4), (2, 3), sharex=ax1s, sharey=ax2s)

ax2.set_yticklabels([])
ax1.set_xticklabels([])
ax2s.set_yticklabels([])
# add minor ticks to the bottom
ax1s.yaxis.set_minor_locator(AutoMinorLocator())
ax2s.yaxis.set_minor_locator(AutoMinorLocator())

### Targeted ###
################
with hdf.File('./buzzard_targetedRealistic_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d', 'ML_pred_3d']
# filter bad values
mask = (target['ML_pred_1d'] != 0)
target = target[mask]

# plot one to one lines
ax1.plot([12, 15.5], [12, 15.5], c='k', zorder=0)
ax2.plot([12, 15.5], [12, 15.5], c='k', zorder=0)
ax3.plot([12, 15.5], [12, 15.5], c='k', zorder=0)
ax4.plot([12, 15.5], [12, 15.5], c='k', zorder=0)
ax1s.axhline(0)
ax2s.axhline(0)
ax3s.axhline(0)
ax4s.axhline(0)

# now for the plotting
###################
#### Power Law ####
###################

for d, c, style, zo in zip([target], ['#7A68A6'], ['-'], [1]):

    print('power law')
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        pyl.log10(d['MASS']),
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax1.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)

    if not c == '#e24a33':
        ax1.fill_between(y_[0],
                         quants[:, 2],
                         quants[:, 0],
                         facecolor=c,
                         alpha=0.3,
                         edgecolor=c)
    err = calc_err(d['MASS'], d['M200c'])
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        err,
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax1s.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)

    if not c == '#e24a33':
        ax1s.fill_between(y_[0],
                          quants[:, 2],
                          quants[:, 0],
                          facecolor=c,
                          alpha=0.3,
                          edgecolor=c)

    print('MAE', median_absolute_error(
        pyl.log10(d['M200c']), pyl.log10(d['MASS'])))
    print('RMSE', pyl.sqrt(mean_squared_error(
        pyl.log10(d['M200c']), pyl.log10(d['MASS']))))

    ############
    #### 1d ####
    ############
    print('1d')
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        d['ML_pred_1d'],
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax2.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax2.fill_between(y_[0],
                         quants[:, 2],
                         quants[:, 0],
                         facecolor=c,
                         alpha=0.4,
                         edgecolor=c)
    err = calc_err(10**d['ML_pred_1d'], d['M200c'])
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        err,
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax2s.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax2s.fill_between(y_[0],
                          quants[:, 2],
                          quants[:, 0],
                          facecolor=c,
                          alpha=0.4,
                          edgecolor=c)

    print('MAE', median_absolute_error(pyl.log10(d['M200c']), d['ML_pred_1d']))
    print('RMSE', pyl.sqrt(mean_squared_error(
        pyl.log10(d['M200c']), d['ML_pred_1d'])))

    #############
    #### 2d #####
    #############
    print('2d')
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        d['ML_pred_2d'],
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax3.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax3.fill_between(y_[0],
                         quants[:, 2],
                         quants[:, 0],
                         facecolor=c,
                         alpha=0.4,
                         edgecolor=c)
    err = calc_err(10**d['ML_pred_2d'], d['M200c'])
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        err,
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax3s.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax3s.fill_between(y_[0],
                          quants[:, 2],
                          quants[:, 0],
                          facecolor=c,
                          alpha=0.4,
                          edgecolor=c)

    print('MAE', median_absolute_error(pyl.log10(d['M200c']), d['ML_pred_2d']))
    print('RMSE', pyl.sqrt(mean_squared_error(
        pyl.log10(d['M200c']), d['ML_pred_2d'])))

    ##############
    ##### 3d #####
    ##############
    print('3d')
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        d['ML_pred_3d'],
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax4.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax4.fill_between(y_[0],
                         quants[:, 2],
                         quants[:, 0],
                         facecolor=c,
                         alpha=0.4,
                         edgecolor=c)
    err = calc_err(10**d['ML_pred_3d'], d['M200c'])
    y_ = astStats.runningStatistic(
        pyl.log10(d['M200c']),
        err,
        pyl.percentile,
        binNumber=20,
        q=[16, 50, 84])
    quants = pyl.array(y_[1])
    ax4s.plot(y_[0], quants[:, 1], style, c=c, zorder=zo)
    if not c == '#e24a33':
        ax4s.fill_between(y_[0],
                          quants[:, 2],
                          quants[:, 0],
                          facecolor=c,
                          alpha=0.4,
                          edgecolor=c)

    print('MAE', median_absolute_error(pyl.log10(d['M200c']), d['ML_pred_3d']))
    print('RMSE', pyl.sqrt(mean_squared_error(
        pyl.log10(d['M200c']), d['ML_pred_3d'])))
    print '----'

### Add Legend ###
##################
line1 = pyl.Line2D([], [], ls='-', color='#7A68A6')
line2 = pyl.Line2D([], [], ls='--', color='#188487')
line3 = pyl.Line2D([], [], ls='-.', color='#e24a33')
#ax1.legend(line1, 'Targeted', loc=2)

#### tweak ####
ax1.set_xticks([12, 13, 14, 15])
ax2.set_xticks([12, 13, 14, 15])
ax2s.set_xticks([12, 13, 14, 15])
ax2s.set_ylim(-2, 4)
ax1s.set_ylim(-2, 4)
ax2s.set_yticks([-2, 0, 2])
ax1s.set_yticks([-2, 0, 2])
ax1.set_ylim(ax2.get_ylim())
ax1s.set_ylim(ax2s.get_ylim())

ax1.set_ylabel('Log $M_{pred}$')
ax1s.set_ylabel('$\epsilon$')
ax1s.set_xlabel('Log $M_{200c}$', fontsize=18)
ax2s.set_xlabel('Log $M_{200c}$', fontsize=18)
ax3s.set_xlabel('Log $M_{200c}$', fontsize=18)
ax4s.set_xlabel('Log $M_{200c}$', fontsize=18)

ax1.text(14, 12.25, 'Power Law', fontsize=18, horizontalalignment='center')
ax2.text(14, 12.25, '$ML_{\sigma}$', fontsize=18, horizontalalignment='center')
ax3.text(14,
         12.25,
         '$ML_{\sigma, z}$',
         fontsize=18,
         horizontalalignment='center')
ax4.text(14,
         12.25,
         '$ML_{\sigma, z, Ngal}$',
         fontsize=18,
         horizontalalignment='center')

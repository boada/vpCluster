import pylab as pyl
from astLib import astStats
import h5py as hdf
from matplotlib.ticker import AutoMinorLocator

def calc_err(pred, true):
    return (pred - true)/true

golden_mean = (pyl.sqrt(5.)-1.0)/2.0
f = pyl.figure(figsize=(7*golden_mean,7))

ax1 = pyl.subplot2grid((3,1), (0,0), rowspan=2)

# now for the bottom bits
ax1s = pyl.subplot2grid((3,1), (2,0))

ax1.set_xticklabels([])
# add minor ticks to the bottom
ax1s.yaxis.set_minor_locator(AutoMinorLocator())

### Targeted ###
################
with hdf.File('./buzzard_targetedRealistic_shifty_masses.hdf5', 'r') as f:
    dset = f[f.keys()[0]]
    target = dset['M200c', 'MASS', 'ML_pred_1d', 'ML_pred_2d2', 'ML_pred_3d']
# filter bad values
mask = (target['ML_pred_1d'] != 0)
target = target[mask]

# plot one to one lines
ax1.plot([11.5,15.5], [11.5,15.5], c='k', zorder=0)
ax1s.axhline(0)

# now for the plotting
###################
#### Power Law ####
###################

c = '#cf4457'
zo = 0

print('power law')
y_ = astStats.runningStatistic(pyl.log10(target['M200c']),
        pyl.log10(target['MASS']), pyl.percentile, binNumber=15, q=[16, 50, 84])
quants = pyl.array(y_[1])
ax1.plot(y_[0],quants[:,1], '-', c=c, zorder=zo)

if not c == '#e24a33':
    ax1.fill_between(y_[0], quants[:,2], quants[:,0], facecolor=c,
        alpha=0.3, edgecolor=c)
err = calc_err(target['MASS'], target['M200c'])
y_ = astStats.runningStatistic(pyl.log10(target['M200c']), err,
        pyl.percentile, binNumber=15, q=[16, 50, 84])
quants = pyl.array(y_[1])
ax1s.plot(y_[0],quants[:,1], '-', c=c, zorder=zo)

if not c == '#e24a33':
    ax1s.fill_between(y_[0], quants[:,2], quants[:,0], facecolor=c,
        alpha=0.3, edgecolor=c)

##############
##### 3d #####
##############

c = '#467821'
zo = 1

print('2d2')
y_ = astStats.runningStatistic(pyl.log10(target['M200c']), target['ML_pred_2d2'],
        pyl.percentile, binNumber=15, q=[16, 50, 84])
quants = pyl.array(y_[1])
ax1.plot(y_[0],quants[:,1], '--', c=c, zorder=zo)
if not c == '#e24a33':
    ax1.fill_between(y_[0], quants[:,2], quants[:,0], facecolor=c,
        alpha=0.4, edgecolor=c)
err = calc_err(10**target['ML_pred_2d2'], target['M200c'])
y_ = astStats.runningStatistic(pyl.log10(target['M200c']), err,
        pyl.percentile, binNumber=15, q=[16, 50, 84])
quants = pyl.array(y_[1])
ax1s.plot(y_[0],quants[:,1], '--', c=c, zorder=zo)
if not c == '#e24a33':
    ax1s.fill_between(y_[0], quants[:,2], quants[:,0], facecolor=c,
        alpha=0.4, edgecolor=c)


### Add Legend ###
##################
line1 = pyl.Line2D([], [], ls='-', color='#cf4457')
line2 = pyl.Line2D([], [], ls='--', color='#467821')
ax1.legend([line1, line2], ['Power Law', '$ML_{\sigma, Ngal}$'], loc=2)

#### tweak ####
ax1.set_xticks([12,13,14,15])
ax1s.set_xticks([12,13,14,15])
ax1s.set_ylim(-2,4)
ax1s.set_yticks([-2,0,2])

ax1.set_ylabel('Log $M_{pred}$')
ax1s.set_ylabel('$\epsilon$')
ax1s.set_xlabel('Log $M_{200c}$', fontsize=18)

#ax1.text(14, 12.25, 'Power Law', fontsize=18, horizontalalignment='center')
#ax4.text(14, 12.25, '$ML_{\sigma, z, Ngal}$', fontsize=18,
#        horizontalalignment='center')
pyl.show()

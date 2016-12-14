import pandas as pd
import pylab as pyl
from glob import glob

files = glob('*.csv')
s = 0
for f in files:
    results = pd.read_csv(f)
    # redshifts
    try:
        q0 = pyl.append(q0, results[results.Q == 0].r.values)
        q1 = pyl.append(q1, results[results.Q == 1].r.values)
        q2 = pyl.append(q2, results[results.Q == 2].r.values)
    except NameError:
        q0 = results[results.Q == 0].r.values
        q1 = results[results.Q == 1].r.values
        q2 = results[results.Q == 2].r.values

# make a figure
f = pyl.figure(1, figsize=(5, 5 * (pyl.sqrt(5.) - 1.0) / 2.0))
ax = f.add_subplot(111)

bins = pyl.linspace(14, 22, 15)
ax.hist(q2,
        weights=pyl.zeros_like(q2) + 1. / q2.size,
        histtype='step',
        bins=bins,
        lw=2,
        label='Q=2')
ax.hist(q1,
        weights=pyl.zeros_like(q1) + 1. / q1.size,
        histtype='step',
        bins=bins,
        lw=2,
        label='Q=1')
ax.hist(q0,
        weights=pyl.zeros_like(q0) + 1. / q0.size,
        histtype='step',
        bins=bins,
        lw=2,
        label='Q=0')

ax.legend(loc='upper right')
ax.invert_xaxis()

ax.set_ylim(0, 0.5)
ax.set_xlabel('$m_r$')
ax.set_ylabel('Fraction of Total')
pyl.show()

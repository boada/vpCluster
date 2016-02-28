import pandas as pd
import pylab as pyl
from glob import glob

files = glob('*.csv')

for f in files:
    results = pd.read_csv(f)
    # good redshifts
    try:
        q0 = pyl.append(q0, results[results.Q == 0].r.values)
        q1 = pyl.append(q1, results[results.Q == 1].r.values)
        x = ~pyl.isnan(results.fiber) & pyl.isnan(results.Q)
        q2 = pyl.append(q2, results.r.values[x.values])
    except NameError:
        q0 = results[results.Q==0].r.values
        q1 = results[results.Q==1].r.values
        x = ~pyl.isnan(results.fiber) & pyl.isnan(results.Q)
        q2 = results.r.values[x.values]

bins = pyl.linspace(14,22,15)
pyl.hist(q2, weights=pyl.zeros_like(q2)+1./q2.size, histtype='step', bins=bins,
        lw=2, label='Q=2')
pyl.hist(q1, weights=pyl.zeros_like(q1)+1./q1.size, histtype='step', bins=bins,
        lw=2, label='Q=1')
q0 = q0[~pyl.isnan(q0)]
pyl.hist(q0, weights=pyl.zeros_like(q0)+1./q0.size, histtype='step', bins=bins,
        lw=2, label='Q=0')

pyl.legend(loc='upper right')
pyl.gca().invert_xaxis()

pyl.ylim(0,0.5)
pyl.xlabel('$m_r$')


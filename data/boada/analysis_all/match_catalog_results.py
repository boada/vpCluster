import numpy as np
import pandas as pd

cluster = 'c203p83+41p0'

results = ['results_may12', 'results_august12', 'results_may13']


cat = pd.read_csv('./catalogs/'+cluster+'_complete.csv')

for result in results:
    print result
    redshifts_part = np.genfromtxt('./redshifts/'+result, names=True, dtype=None)

    # filter out for the cluster that we want
    mask = redshifts_part['cluster'] == cluster
    redshifts_part = redshifts_part[mask]
    redshifts_part = pd.DataFrame(redshifts_part)

    try:
        redshifts = redshifts.append(redshifts_part)
    except NameError:
        redshifts = redshifts_part

merged = pd.merge(cat, redshifts, left_on=['tile', 'dither', 'fiber'],
        right_on=['tile','dither','fiber'], how='outer')

del merged['index']

merged.to_csv(cluster+'redshifts.csv', index=False)


import pandas as pd
from glob import glob
import numpy as np
from itertools import combinations

# get the files
files = glob('../redshifts/*_redshifts.csv')

# load the frames into an array
frames = [pd.read_csv(f) for f in files]

# concat all of the frames together
result = pd.concat(frames, ignore_index=True)

# find the results that aren't nan
mask = ~np.isnan(result.Q)
result = result[mask]

# build the features
features = pd.DataFrame()
# mags first
features = result[['u','g','r', 'i', 'z']]

# colors
colors = combinations('zirgu', 2)
for i in colors:
    features['%s-%s' % (i[0], i[1])] = result.loc[:,i[0]] - result.loc[:,i[1]]

# colors squared
colors = combinations('zirgu', 2)
for i in colors:
    features['(%s-%s)2' % (i[0], i[1])] = (result.loc[:,i[0]] -
            result.loc[:,i[1]])**2






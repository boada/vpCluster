import scipy.stats as stats
from glob import glob
import pandas as pd
from numpy import sort as sort

#files = glob('*members.csv')
#for f in files:
def mod_dom(f):
    print f
    data = pd.read_csv(f)
    f = f.split('/')[1]
    data2 = pd.read_csv('catalogs/'+ f.split('_')[0] + '_complete.csv')

    c = pd.merge(data, data2, left_on='ra', right_on='ra', how='inner')

    kur = stats.kurtosis(c['redshift'])
    ske = stats.skew(c['redshift'])
    mod = (1 + ske**2)/(3 + kur**2)

    mags = c['r'].values
    mags = sort(mags)
    dom = abs(mags[0] - mags[1])

    #print mod, dom
    return mod, dom




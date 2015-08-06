import pandas as pd
from glob import glob
import numpy as np
from astLib import astStats as ast
from astLib import astCalc as aca
import sys

def findClusterCenterRedshift(data, errors=False):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator. Puts the result in the ClusZ column of the data
    array. If errors = True, then this will return the 95% C.I. from 1000
    bootstrap shuffles.

    '''

    x = np.copy(data.redshift[data.interloper == 'NO'].values)
    avgz = ast.biweightLocation(x, tuningConstant=6.0)
    #return ast.biweightClipped(data['Z'], 6.0, 3)['biweightLocation']
    #return ast.biweightLocation(data['Z'], tuningConstant=6.0)

    print len(x)
    if errors:
        ci = ast.bootstrap(x, ast.biweightLocation, tuningConstant=6.0)
        return avgz, ci
    else:
        return avgz

def calcLOSVD(data, errors=False):
    ''' Using the previously computed LOSVs we will find the LOSVD. This will
    give us a few options to do that based on the number of objects that we
    have. If errors = True, then this will return the 95% C.I. from 1000
    bootstrap shuffles.

    '''

    if data.interloper.value_counts().NO >= 15:
        x = np.copy(data.LOSV[data.interloper=='NO'].values)
        LOSVD = ast.biweightScale_test(x, tuningConstant=9.0)
        if errors:
            ci = ast.bootstrap(x, ast.biweightScale_test, tuningConstant=9.0)
            return LOSVD, ci
        else:
            return LOSVD

    else:
        x = np.copy(data.LOSV[data.interloper=='NO'].values)
        LOSVD = ast.gapperEstimator(x)
        if errors:
            ci = ast.bootstrap(x, ast.gapperEstimator)
            return LOSVD, ci
        else:
            return LOSVD

def calc_mass_Evrard(data, A1D=1177., alpha=0.364, errors=False):
    ''' This uses the relation from Munari2013 to calculate the halo mass from
    the observed velocity dispersion. The chosen scaling relations are from
    their table 1 which has been calibrated using galaxies and not dark matter
    halos only.

    '''

    avgz = findClusterCenterRedshift(data)
    vd = calcLOSVD(data)

    if errors:
        mass = 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)
        ci = ast.bootstrap(data, calc_mass_Evrard)
        return mass, ci
    else:
        return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        clusters = glob('./members/*_members.csv')
        #clusters = [c.rstrip('_members.csv') for c in clusters]

        for c in clusters:
            data = pd.read_csv(c)
            mass = calc_mass_Evrard(data)
            c = c.rstrip('_members.csv')
            print c, mass/1e15

    elif len(sys.argv) == 2:
        cluster = './members/'+sys.argv[1]+'_members.csv'
        data = pd.read_csv(cluster)
        mass = calc_mass_Evrard(data)
        cluster = cluster.rstrip('_members.csv')
        print cluster, mass/1e15

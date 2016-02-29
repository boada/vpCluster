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
    w = np.copy(data.redshift_err.values) * 2.
    #avgz = ast.biweightLocation(x, tuningConstant=6.0)
    #return ast.biweightClipped(data['Z'], 6.0, 3)['biweightLocation']
    #return ast.biweightLocation(data['Z'], tuningConstant=6.0)
    avgz = np.average(x, weights=1./w)

    #print len(x)
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

def calc_mass_Evrard(data, A1D=1177., alpha=0.364):
    ''' This uses the relation from Munari2013 to calculate the halo mass from
    the observed velocity dispersion. The chosen scaling relations are from
    their table 1 which has been calibrated using galaxies and not dark matter
    halos only.

    '''

    avgz = findClusterCenterRedshift(data)
    vd = calcLOSVD(data)

    if avgz == None:
        pass
    else:
        return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)
        #return 1/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

if __name__ == "__main__":
    ''' This takes the membership information from the member catalogs and
    computes the LOSVD and cluster mass from the member galaxies. Right now it
    is using a simple power law to do this, but things will change in the
    future when I add the stuff from the des study.

    '''


    if len(sys.argv) == 1:
        clusters = glob('./members/*_members.csv')
        #clusters = [c.rstrip('_members.csv') for c in clusters]

        alpha = 0.32 # 95% CI.
        resamples = 500
        with open('cluster_props', 'w') as f:
            f.write('#name mass mass_lower mass_upper\n')
            for c in clusters:
                data = pd.read_csv(c)
                data = data[data.interloper=='NO']
                mass = calc_mass_Evrard(data)

                # this reasamples the dataframe to bootstrap the error on the mass
                idx = np.random.random_integers(0, len(data)-1,
                        size=(resamples, len(data)-1))
                mass_ci = np.sort([calc_mass_Evrard(data.iloc[row]) for row in
                    idx])
                limits = (mass_ci[int((alpha/2.0) * resamples)],
                        mass_ci[int((1-alpha/2.0) * resamples)])

                c = c.rstrip('_members.csv')
                try:
                    print c, mass, limits[0], limits[1]
                except TypeError:
                    print c, mass/1e15, '---', limits[1]/1e15

                f.write(c.split('/')[-1]+' '+str(mass)+' '+str(limits[0])+\
                    ' '+str(limits[1])+'\n')


    elif len(sys.argv) == 2:
        cluster = './members/'+sys.argv[1]+'_members.csv'
        data = pd.read_csv(cluster)
        mass = calc_mass_Evrard(data)
        cluster = cluster.rstrip('_members.csv')
        print cluster, mass/1e15

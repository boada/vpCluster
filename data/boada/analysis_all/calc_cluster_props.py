import pandas as pd
import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
import sys

def findClusterCenterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator. Puts the result in the ClusZ column of the data
    array.

    '''

    data['CLUSZ'] = ast.biweightLocation(data['Z'], tuningConstant=6.0)
    #return ast.biweightClipped(data['Z'], 6.0, 3)['biweightLocation']
    #return ast.biweightLocation(data['Z'], tuningConstant=6.0)

    return data

def calc_mass_Evrard(data, A1D=1177., alpha=0.364):
    ''' This uses the relation from Munari2013 to calculate the halo mass from
    the observed velocity dispersion. The chosen scaling relations are from
    their table 1 which has been calibrated using galaxies and not dark matter
    halos only.

    '''

    avgz = data.Clusz
    vd = data.LOSV

    return 1e15/(aca.H0 * aca.Ez(avgz)/100.) * (vd/A1D)**(1/alpha)

def main(cluster):
    pass

if __name__ == "__main__":
    main(sys.argv[1])

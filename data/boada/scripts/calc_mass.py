from math import sqrt
from astLib import astStats
from astLib import astCalc
import pandas as pd
from glob import glob

def calcVD_big(data):

    return astStats.biweightScale(data, tuningConstant=9.0)

def calcVD_small(data):

    return astStats.gapperEstimator(data)

def calc_mass(data):

    if len(data) > 10:
        vd = calcVD_big
    else:
        vd = calcVD_small

    avgz = astStats.biweightLocation(data, tuningConstant=6.0)

    r200 = sqrt(3) * vd /(10*astCalc.H0 * astCalc.Ez(avgz))

    m200 = 3 * vd**2 * r200 * 1000**3 * 3.08E19/6.67384E-11

    return m200

if __name__ == "__main__":
    files = glob('members.csv')
    for f in files
        data = pd.read_csv(f)
        mass = calc_mass(data['LOSV'].values)
        print len(data), mass/1.9891E30

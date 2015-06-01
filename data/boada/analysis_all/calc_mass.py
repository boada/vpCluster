from math import sqrt
from astLib import astStats
from astLib import astCalc
import pandas as pd
from glob import glob
from bootstrap import bootstrap
def calcVD_big(data):

    return astStats.biweightScale(data, tuningConstant=9.0)

def calcVD_small(data):

    return astStats.gapperEstimator(data)

def calc_mass(data):

    if len(data) > 10:
        vd = calcVD_big(data['LOSV'].values)
        up, low = bootstrap(data['LOSV'].values, astStats.biweightScale,
                alpha=0.32, tuningConstant=9.0)
    else:
        vd = calcVD_small(data['LOSV'].values)
        up, low = bootstrap(data['LOSV'].values, astStats.gapperEstimator,
                alpha=0.32)

#    print vd, abs(vd-up), abs(vd-low),

    avgz = astStats.biweightLocation(data['redshift'].values, tuningConstant=6.0)

    r200 = sqrt(3) * vd /(10*astCalc.H0 * astCalc.Ez(avgz))
    r200up = sqrt(3) * up /(10*astCalc.H0 * astCalc.Ez(avgz))
    r200low = sqrt(3) * low /(10*astCalc.H0 * astCalc.Ez(avgz))

#    print r200, abs(r200-r200up), abs(r200-r200low),

    a = 3 * sqrt(3) * 1000**3 * 3.08E19/(10*astCalc.H0 * astCalc.Ez(avgz) *\
            6.67384E-11)

    m200 = a * vd**3

    # propagate errors
    m200low = a * vd**2 * 3 * low
    m200up = a * vd**2 * 3 * up

    print vd, m200/1.9891E30, m200low/1.9891E30, m200up/1.9891E30

    return data, m200, vd

files = glob('*members.csv')
for f in files:
    data = pd.read_csv(f)
    print f.split('_')[0],
    mass = calc_mass(data)
    #print len(data), mass/1.9891E30

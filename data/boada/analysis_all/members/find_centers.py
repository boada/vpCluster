import pandas as pd
import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
import sys

c = 2.99E5  # speed of light in km/s


def findClusterCenterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator.

    '''

    # this filters out interlopers and copies in one step
    x = np.copy(data.redshift[data.interloper == 'NO'].values)
    #x = np.copy(data['redshift'].values)
    return ast.biweightLocation(x, tuningConstant=6.0)


def findseparationSpatial(data, center):
    ''' Finds the distance to all of the galaxies from the center of the
    cluster in the spatial plane. Returns values in Mpc.

    '''

    # Add a new column to the dataframe
    data['separation'] = 0.0
    for row in data.iterrows():
        sepDeg = aco.calcAngSepDeg(center[0], center[1], row[1]['ra'],
                                   row[1]['dec'])
        sepMpc = sepDeg * aca.da(row[1]['redshift']) / 57.2957795131
        data['separation'][row[0]] = sepMpc
        #data['separation'][row[0]] = sepDeg * 3600

    return data


def findLOSV(data, avgz):
    ''' Finds the line of sight velocity for each of the galaxies.

    '''

    c = 2.99E5  # speed of light in km/s

    # Add a new column to the dataframe
    data['LOSV'] = 0.0
    for row in data.iterrows():
        data['LOSV'][row[0]] = c * (row[1]['redshift'] - avgz) / (1 + avgz)

    return data


def split_list(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list will have an extra element.

    '''
    length = len(alist)
    return [alist[i * length // wanted_parts:(i + 1) * length // wanted_parts]
            for i in range(wanted_parts)]


def rejectInterlopers(data):
    ''' Does all of the work to figure out which galaxies don't belong. Makes
    several sorted copies of the dataframe and then applies the fixed gapper
    method.

    '''

    # mask out the interlopers
    data_orig = data
    data = data[data.interloper == 'NO']

    # make some copies so we can sort them around
    sepSorted = data.sort('separation', ascending=True)
    # How many parts to break into
    parts = len(data) // 10
    splitData = split_list(sepSorted, parts)

    # Now we sort the parts by LOSV and find the rejects
    interlopers = []
    for part in splitData:
        # sort by LOSV
        LOSVsorted = part.sort('LOSV', ascending=True)
        rejected = True
        while rejected:
            # Find the difference between all of the neighboring elements
            difference = np.diff(LOSVsorted['LOSV'])
            # If diff > 1000 reject
            rejects = abs(difference) > 1000
            # Now remove those items
            indices = np.where(rejects == True)
            #print LOSVsorted['LOSV']
            #print difference
            #print indices[0]
            if rejects.any() == True:
                # Always take the more extreme index
                for index, i in enumerate(indices[0]):
                    if (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                            abs(LOSVsorted['LOSV'][LOSVsorted.index[i + 1]])
                        ) > 0:
                        pass
                    elif (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                          abs(LOSVsorted['LOSV'][LOSVsorted.index[i + 1]])
                          ) < 0:
                        indices[0][index] = i + 1

                #print LOSVsorted.index[list(indices[0])]
                dataframeIndex = list(LOSVsorted.index[list(indices[0])])
                LOSVsorted = LOSVsorted.drop(dataframeIndex)
                interlopers += dataframeIndex
            else:
                rejected = False
    print 'interlopers', interlopers

    data['interloper'][interlopers] = 'YES'

    #return data.drop(interlopers)
    #return data_orig.update(data)
    data_orig.update(data)


def rejectInterlopers_group(data, avgz, sigmav=500):

    deltaZmax = 2 * sigmav * (1. + avgz) / c
    deltaRmax = (c * deltaZmax) / (9.5 * aca.H0 * aca.Ez(avgz))  # Mpc
    #deltaThetamax = 206265 * deltaRmax * aca.da(avgz) #arcseconds

    # make redshift cut
    data.interloper[abs(data.redshift - avgz) >= deltaZmax] = 'YES'
    #data = data[abs(data.redshift - avgz) <= deltaZmax]

    # make spatial cut
    data.interloper[data.separation >= deltaRmax] = 'YES'
    #data = data[data.separation <= deltaRmax]

    return data


if __name__ == "__main__":
    ''' I don't really remember what this script is supposed to do but it looks
    like it it has something to do with the miscentering problem that we were
    trying to look at originally. If I had a guess it computes the center of
    the cluster based on the member galaxies and them compares it to the center
    that we have given in the centers list.

    This script isn't to be used for anything at the momement.

    '''

    catalog = './../redshifts/' + sys.argv[1] + '_redshifts.csv'
    # get the center
    with open('./../centers/' + sys.argv[1] + '_center.list', 'r') as f:
        line = f.read()
    line = line.split(' ')
    center = float(line[0]), float(line[1])

    results = pd.read_csv(catalog)
    # mask out the missing values
    mask = ~np.isnan(results['redshift'])
    results = results[mask]

    # now we start the process
    separated = findseparationSpatial(results, center)
    # add column for interlopers

    galaxyNum = separated.shape[0]

    if galaxyNum >= 20:
        while True:
            try:
                avgz = findClusterCenterRedshift(cleaned)
                losv = findLOSV(cleaned, avgz)
                losv['interloper'] = 'NO'
                losv = findseparationSpatial(losv, (RAcenter, DECcenter))
            except NameError:
                separated['interloper'] = 'NO'
                members = separated.interloper.value_counts().NO
                avgz = findClusterCenterRedshift(separated)
                #avgz = 0.2256
                losv = findLOSV(separated, avgz)

            rejectInterlopers(losv)
            cleaned = losv

            print avgz
            if members == cleaned.interloper.value_counts().NO:
                break
            else:
                members = cleaned.interloper.value_counts().NO

            mems = cleaned[cleaned.interloper == 'NO']
            #Mr = [aca.absMag(r, z) for r, z in zip(mems.r, mems.redshift)]

            RAcenter = np.average(mems.ra, weights=1. / np.array(mems.r))
            DECcenter = np.average(mems.dec, weights=1. / np.array(mems.r))
    else:
        while True:
            try:
                avgz = findClusterCenterRedshift(cleaned)
                losv = findLOSV(cleaned, avgz)
                losv['interloper'] = 'NO'
            except NameError:
                separated['interloper'] = 'NO'
                members = separated.interloper.value_counts().NO
                avgz = findClusterCenterRedshift(separated)
                losv = findLOSV(separated, avgz)

            try:
                rejectInterlopers_group(losv, avgz, sigmav=LOSVD)
                cleaned = losv
                LOSVD = 1.135 *\
                ast.gapperEstimator(cleaned.LOSV[cleaned.interloper == 'NO'])
                print LOSVD
            except NameError:
                rejectInterlopers_group(losv, avgz)
                cleaned = losv
                LOSVD = 1.135 *\
                ast.gapperEstimator(cleaned.LOSV[cleaned.interloper == 'NO'])
                print LOSVD

            print avgz
            if members == cleaned.interloper.value_counts().NO:
                break
            else:
                members = cleaned.interloper.value_counts().NO

    mems = cleaned[cleaned.interloper == 'NO']
    #Mr = [aca.absMag(r, z) for r, z in zip(mems.r, mems.redshift)]

    RAcenter = np.average(mems.ra, weights=1. / np.array(mems.r))
    DECcenter = np.average(mems.dec, weights=1. / np.array(mems.r))

    print center, RAcenter, DECcenter

    #cleaned.to_csv(sys.argv[1]+'_members.csv', index=False)

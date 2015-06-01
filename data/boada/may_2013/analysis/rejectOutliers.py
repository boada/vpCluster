import glob
import pandas as pd
import pylab as pyl
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca

def parseResults(files):
    ''' Reads all of the results files and puts them into a list with the
    results. Returns field, dither, fiber, and redshift.

    '''

    r = []
    for f in files:
        print f
        cluster, field, dither = f.split('_')
        data = pyl.genfromtxt(f, delimiter='\t', names=True, dtype=None)
        try:
            for fiber, z, Q in zip(data['Fiber'], data['Redshift'],
                    data['Quality']):
                if Q == 0:
                    r.append((field, 'D'+str(dither.rstrip('.results')), fiber,
                        z))
        except TypeError:
            fiber = int(data['Fiber'])
            z = float(data['Redshift'])
            Q = int(data['Quality'])
            if Q == 0:
                r.append((field, 'D'+str(dither.rstrip('.results')), fiber, z))

    print len(r), 'objects read'
    return r

def matchToCatalog(results, catalog):
    ''' Matches the list pumped out from parseResults to the full dataframe.
    Returns another pandas dataframe to work with.

    '''

    cat = pd.read_csv(catalog)[['tile', 'dither', 'fiber', 'ra', 'dec']]
    data = pd.DataFrame(results, columns=['tile', 'dither', 'fiber',
        'redshift'])

    # This does the actual matching
    matched = pd.merge(cat, data, left_on=['tile', 'dither', 'fiber'],
            right_on=['tile', 'dither', 'fiber'], how='inner')

    return matched

def findClusterCenterRedshift(data):
    ''' Finds the center of the cluster in redshift space using the
    biweightlocation estimator.

    '''
    x = pyl.copy(data['redshift'].values)
    return ast.biweightLocation(x, tuningConstant=6.0)

def findSeperationSpatial(data, center):
    ''' Finds the distance to all of the galaxies from the center of the
    cluster in the spatial plane. Returns values in Mpc.

    '''

    # Add a new column to the dataframe
    data['seperation'] = 0.0
    for row in data.iterrows():
        sepDeg = aco.calcAngSepDeg(center[0], center[1], row[1]['ra'],
                row[1]['dec'])
        sepMpc = sepDeg * aca.da(row[1]['redshift'])/57.2957795131
        data['seperation'][row[0]] = sepMpc

    return data

def findLOSV(data):
    ''' Finds the line of sight velocity for each of the galaxies.

    '''

    c = 2.99E5 # speed of light in km/s

    avgz = findClusterCenterRedshift(data)

    # Add a new column to the dataframe
    data['LOSV'] = 0.0
    for row in data.iterrows():
        data['LOSV'][row[0]] = c *(row[1]['redshift'] - avgz)/(1 + avgz)

    return data

def split_list(alist, wanted_parts=1):
    ''' Breaks a list into a number of parts. If it does not divide evenly then
    the last list will have an extra element.

    '''
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts]
        for i in range(wanted_parts) ]

def rejectInterlopers(data):
    ''' Does all of the work to figure out which galaxies don't belong. Makes
    several sorted copies of the dataframe and then applies the fixed gapper
    method.

    '''

    # make some copies so we can sort them around
    sepSorted = data.sort('seperation', ascending=True)
    # How many parts to break into
    parts = len(data)//15
    splitData = split_list(data, parts)

    # Now we sort the parts by LOSV and find the rejects
    interlopers = []
    for part in splitData:
        # sort by LOSV
        LOSVsorted = part.sort('LOSV', ascending=True)
        rejected = True
        while rejected:
            # Find the difference between all of the neighboring elements
            difference = pyl.diff(LOSVsorted['LOSV'])
            # If diff > 1000 reject
            rejects = abs(difference) > 1000
            # Now remove those items
            indices = pyl.where(rejects == True)
            #print LOSVsorted['LOSV']
            #print difference
            #print indices[0]
            if rejects.any() == True:
                # Always take the more extreme index
                for index, i in enumerate(indices[0]):
                    if (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                        abs(LOSVsorted['LOSV'][LOSVsorted.index[i+1]])) > 0:
                            pass
                    elif (abs(LOSVsorted['LOSV'][LOSVsorted.index[i]]) -
                        abs(LOSVsorted['LOSV'][LOSVsorted.index[i+1]])) < 0:
                            indices[0][index] = i+1

                #print LOSVsorted.index[list(indices[0])]
                dataframeIndex = list(LOSVsorted.index[list(indices[0])])
                LOSVsorted = LOSVsorted.drop(dataframeIndex)
                interlopers += dataframeIndex
            else:
                rejected = False
    print 'interlopers',interlopers
    return data.drop(interlopers)

catalog = '/Users/steven/Projects/cluster/data/boada/may_2012/catalogs/c260p61+32p13_complete.csv'

files = glob.glob('*.results')
center = 260.61326, 32.132568

results = parseResults(files)
matched = matchToCatalog(results, catalog)
seperated = findSeperationSpatial(matched, center)
losv = findLOSV(seperated)

cleaned = rejectInterlopers(losv)


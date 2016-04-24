import pandas as pd
import numpy as np
from astLib import astCoords as aco
from astLib import astStats as ast
from astLib import astCalc as aca
import sys
from glob import glob

# buzzard simulation cosmology
aca.H0 = 70
aca.OMEGA_M0 = 0.286
aca.OMEGA_L0 = 0.714


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
    sepDeg = np.array(aco.calcAngSepDeg(center[0], center[1], data.ra.values,
            data.dec.values))
    da = np.array([aca.da(z)/57.2957795131 for z in data.redshift.values])

    sepMpc = da * sepDeg
    data.loc[:, 'separation'] = sepMpc

    return data

def findLOSV(data, avgz):
    ''' Finds the line of sight velocity for each of the galaxies. Computes the
    LOSVD error as 2* redshift error.

    '''

    c = 2.99E5 # speed of light in km/s

    LOSV = c * (data.redshift.values - avgz)/ (1+avgz)
    LOSV_err = c/(1+avgz) * data.redshift_err.values * 2

    # Add a new column to the dataframe
    data.loc[:, 'LOSV'] = LOSV
    data.loc[:, 'LOSV_err'] = LOSV_err

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

    # mask out the interlopers
    data_orig = data
    data = data[data.interloper == 'NO']

    # make some copies so we can sort them around
    sepSorted = data.sort_values(by='separation', ascending=True)
    # How many parts to break into
    parts = len(data)//10
    splitData = split_list(sepSorted, parts)

    # Now we sort the parts by LOSV and find the rejects
    interlopers = []
    for part in splitData:
        # sort by LOSV
        LOSVsorted = part.sort_values(by='LOSV', ascending=True)
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

    data.loc[interlopers,'interloper'] = 'YES'

    #return data.drop(interlopers)
    #return data_orig.update(data)
    data_orig.update(data)

def shifty_gapper(r, z, zc, vlimit=10000, ngap=30, glimit=1000):
    '''Determine cluster membersip according to the shifty
      gapper technique
      The shifty gapper technique includes the following steps:
      1.  Remove outliers outside of +- vlimit
      2.  Locate the Ngap targets with radius close to r_i
      3.  Within this sample of N targets, identify sources with
          |v_pec| < |v_pec_i|
      4.  Meaure the largest gap within this subset of targets
      5.  If v_gap is larger than glimit, reject the source
      Parameters
      -----------
      vlimit: float
         Initial limit in velocity offset from the cluster center
      ngap: int
         Number of sources to use in in estimating the gap
      glimit:  float
         Maximum gap size for rejecting objects

      Parameters
      -----------
      incluster: ndarray
          Returns a boolean array where objects in the cluster have
          a value of True
    '''

    def z2v(z, zc):
        '''Convert the redshift to km/s relative to the cluster center '''
        return 2.99792458e5*(z-zc)/(1+zc)

    #convert to the velocity scale
    v = z2v(z,zc)

    #limit the source to only sources within the vlimit
    vmask = abs(v) < vlimit

    nobj=len(r)
    incluster=np.zeros(nobj, dtype=bool)

    if nobj<ngap:
        raise Exception('Number of sources is less thant number of gap sources')

    for i in range(nobj):
        if abs(v[i])<vlimit:
            #find the ngap closest sources
            r_j=abs(r[vmask]-r[i]).argsort()
            vg=v[vmask][r_j[0:ngap]]

            #find the sources with |v_pec| < |v_pec_i|
            mask=abs(vg)<=abs(v[i])
            if mask.sum()>1:
                vg=vg[mask]
                #now sort these sources and find the biggest gap
                vg.sort()
                if np.diff(vg).max()<glimit: incluster[i]=True
            else:
                incluster[i]=True

    return incluster


def rejectInterlopers_group(data, avgz, sigmav=500):

    c = 2.99E5 # speed of light in km/s
    deltaZmax = 2 * sigmav * (1.+avgz) / c
    deltaRmax = (c * deltaZmax)/(9.5*aca.H0*aca.Ez(avgz)) # Mpc
    #deltaThetamax = 206265 * deltaRmax * aca.da(avgz) #arcseconds

    # make redshift cut
    data.loc[abs(data.redshift - avgz) >= deltaZmax, 'interloper'] = 'YES'
    #data.interloper[abs(data.redshift - avgz) >= deltaZmax] = 'YES'
    #data = data[abs(data.redshift - avgz) <= deltaZmax]

    # make spatial cut
    data.loc[data.separation >=deltaRmax, 'interloper'] = 'YES'
    #data.interloper[data.separation >=deltaRmax] = 'YES'
    #data = data[data.separation <= deltaRmax]

    return data

if __name__ == "__main__":
    ''' This script does the interlop rejection of the redshift catalogs.
    Currently, it uses my own version of the shiffting gapper method, but I
    might try adding a few other's implimentations of it to see how things
    might change. This creates new files with the membership information for
    use in the science and for making the tables with membership.

    '''

    if len(sys.argv) == 1:
        clusters = glob('../redshifts/*_redshifts.csv')
        # strip off the extra wording
        clusters = [c.rstrip('_redshifts.csv') for c in clusters]
        clusters = [c.split('/')[2] for c in clusters]

    elif len(sys.argv) == 2:
        clusters = [sys.argv[1]]

    for c in clusters:
        print c
        catalog = './../redshifts/'+c+'_redshifts.csv'
        print catalog
        # get the center
        with open('./../centers/'+c+'_center.list', 'r') as f:
            line = f.read()
        line = line.split(' ')
        center = float(line[0]), float(line[1])
        avgz = float(line[2])

        results = pd.read_csv(catalog)

        # mask out the bad redshifts
        mask = (results.Q != 2)
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
                except NameError:
                    separated.loc[:, 'interloper'] = 'NO'
                    members = separated.interloper.value_counts().NO
                    #avgz = findClusterCenterRedshift(separated)
                    #avgz = 0.2256
                    losv = findLOSV(separated, avgz)

                mask = shifty_gapper(losv.separation.values,
                        losv.redshift.values, avgz, ngap=20, vlimit=5000)
                #rejectInterlopers(losv)
                losv.loc[~mask, 'interloper'] = 'YES'
                cleaned = losv

                print avgz
                if members == cleaned.interloper.value_counts().NO:
                    break
                else:
                    members = cleaned.interloper.value_counts().NO

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

#        cleaned.to_csv(c+'_members.csv', index=False)
        del cleaned

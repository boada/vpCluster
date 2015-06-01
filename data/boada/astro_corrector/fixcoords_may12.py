import pyfits as pyf
from astLib import astSED
from astLib import astWCS
from astLib import astStats
from astLib import astCoords
import numpy as np

def loadIFUSpectra(objectFileName):
    """ Loads in an object spectrum - this has to be in DEEP2 pipeline
    spec1d format (i.e. fits tables)
    Object spectrum is smoothed by boxcar of size smoothPix.

    Returns a dictionary containing object and sky astSED.SED objects
    {'object', 'sky'}

    """

    print "Loading IFU spectrum ..."

    oimg = pyf.open(objectFileName)

    # Load the IFU data -- Row-stacked spectra
    odata = oimg[1].data
    odata_dim = odata.shape
    wcs = astWCS.WCS(objectFileName, extensionName=1)
    owavelengthStartEnd = wcs.getImageMinMaxWCSCoords()[0:2]
    fiberNumber = wcs.getImageMinMaxWCSCoords()[2:4]
    owavelengthStep = oimg[1].header['CDELT1']

    owavelengthRange = [owavelengthStartEnd[0] + i * owavelengthStep
                        for i in range(odata_dim[1])]

    # Check to make sure we got it right
    if not owavelengthRange[-1] == owavelengthStartEnd[-1]:
        print 'The ending wavelenghts do not match... Exiting'
        sys.exit(1)
    else:
        sums = [sum(odata[i,:]) for i in range(odata.shape[0])]
        #find the median value of all the fibers
        med = astStats.clippedMedianStdev(sums)
        med = med['clippedMedian']

        skyfibers = [i for i in range(odata.shape[0])\
                if sum(odata[i,:]) <= med]
        skydata = odata.take(skyfibers, axis=0)

        oskyflux = [np.average(skydata[:,i])\
                for i in range(skydata.shape[1])]

    RSS = []
    for i in range(int(fiberNumber[1])):
        oflux = odata[i] - oskyflux

        # Mask out extreme values in spectrum
        # Just because edges dodgy in efosc
        med = np.median(oflux)
        oflux[np.greater(abs(oflux), 10.0*med)] = 0.0001

        objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)

        #  make it > 0 everywhere
        #objSED.flux = objSED.flux - objSED.flux.min()
        #objSED.flux = objSED.flux / objSED.flux.max()

        skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)

        RSS.append({'object': objSED, 'sky': skySED})
    return RSS

def find_center(points, weights, star):
    x = [p[0] for p in points]
    y = [p[1] for p in points]

    c = (np.average(x, weights=weights), np.average(y, weights=weights))
    #print w
    print c

    Ra_shift = astCoords.calcAngSepDeg(c[0],c[1],star[0],c[1])*3600
    Dec_shift = astCoords.calcAngSepDeg(c[0],c[1],c[0],star[1])*3600

# Need to shift to the South West, so make the shifts negative
    new_coords = astCoords.shiftRADec(c[0], c[1], -1*Ra_shift, -1*Dec_shift)
    print Ra_shift, Dec_shift
    print new_coords

def mk_weights(fibers, RSS):
    w = []
    for f in fibers:
        d = [flux for flux, wl in zip(RSS[246-f]['object'].flux,
            RSS[246-f]['object'].wavelength) if 5000<wl<5010]
        w.append(sum(d))

    return w

def load_coords(fibers, fiberPositions):
    # determines which fibers should be included and returns the fiber number
    # along with the position (RA,DEC) of each.

    coords = np.genfromtxt(fiberPositions, dtype=str)

    temp = []
    for f in fibers:
        ra = coords[f-1, 1]
        dec = coords[f-1,2]
        ra_deg = astCoords.hms2decimal(ra, ':')
        dec_deg = astCoords.dms2decimal(dec, ':')
        #temp.append({'fiber':f, 'ra':ra_deg, 'dec':dec_deg})
        temp.append((ra_deg, dec_deg))

    return temp


#star = (203.861012858877,41.0105477333859)
#star = (203.84062, 41.01396)
#star = (203.80628, 41.02194)
star = (203.81698, 40.98942)

#fibersD1 = 129, 130, 143, 144, 145, 158, 159
#fibersD2 = 129, 130, 143, 144, 145, 158, 159
#fibersD3 = 129, 130, 143, 144, 145, 158, 159
#fibersD1 = fibersD2 = fibersD3 = 93,107,108,122

#fibersD1 = fibersD2 = fibersD3 = 35,36,49,50

fibersD1 = fibersD2 = fibersD3 = 82,83,96,97,98,111,112

D1= './data/bcs0048_crcl_oextr1.fits'
D2= './data/bcs0049_crcl_oextr1.fits'
D3= './data/bcs0184_crcl_oextr1.fits'

pointsD1 = load_coords(fibersD1,'./coords/C203p83+41p00_SW_D1_coords.txt')
pointsD2 = load_coords(fibersD2,'./coords/C203p83+41p00_SW_D2_coords.txt')
pointsD3 = load_coords(fibersD3,'./coords/C203p83+41p00_SW_D3_coords.txt')

RSSD1 = loadIFUSpectra(D1)
RSSD2 = loadIFUSpectra(D2)
RSSD3 = loadIFUSpectra(D3)

wD1 = mk_weights(fibersD1, RSSD1)
wD2 = mk_weights(fibersD2, RSSD2)
wD3 = mk_weights(fibersD3, RSSD3)

center = find_center(pointsD1 + pointsD2 + pointsD3, wD1+wD2+wD3, star)



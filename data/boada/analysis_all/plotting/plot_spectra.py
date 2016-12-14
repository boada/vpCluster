import pyfits
from astLib import astWCS, astSED
import pylab as pyl
import sys


def loadFile(objectFileName):
    oimg = pyfits.open(objectFileName)

    # Load the IFU data -- Row-stacked spectra
    odata = oimg[1].data
    oError = oimg[2].data
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
        # make median sky
        specs = pyl.array([flux for flux in odata])
        skySpec = pyl.median(specs, axis=0)

    RSS = []
    for i in range(int(fiberNumber[1])):
        #oflux = odata[i] - oskyflux
        oflux = odata[i] - skySpec
        oflux[pyl.isnan(oflux)] = 0.0
        oErrorFlux = oError[i]
        #oflux = odata[i]

        # Mask out extreme values in spectrum
        # Just because edges dodgy in efosc
        med = pyl.median(oflux)
        oflux[pyl.greater(abs(oflux), 10.0 * med)] = 0.0001

        objSED = astSED.SED(wavelength=owavelengthRange, flux=oflux)
        #skySED = astSED.SED(wavelength=owavelengthRange, flux=oskyflux)
        skySED = astSED.SED(wavelength=owavelengthRange, flux=skySpec)
        errSED = astSED.SED(wavelength=owavelengthRange, flux=oErrorFlux)

        #  make it > 0 everywhere
        objSED.flux = objSED.flux - objSED.flux.min()
        objSED.flux = objSED.flux / objSED.flux.max()
        errSED.flux = errSED.flux - errSED.flux.min()
        errSED.flux = errSED.flux / errSED.flux.max()
        skySED.flux = skySED.flux - skySED.flux.min()
        skySED.flux = skySED.flux / skySED.flux.max()

        RSS.append({'object': objSED, 'sky': skySED, 'error': errSED})

    return RSS


filepath = './../../may_2012/analysis/c203p83+41p0/'

specs = [{'file': 'c203p83+41p0_NE_3.fits',
          'fiber': 106,
          'z': 0.22354}, {'file': 'c203p83+41p0_SE_3.fits',
                          'fiber': 198,
                          'z': 0.11766}, {'file': 'c203p83+41p0_NW_3.fits',
                                          'fiber': 20,
                                          'z': 0.2310}]

# make some figures
f, ax = pyl.subplots(3, 1, squeeze=True)

for i, (spec, c) in enumerate(zip(specs, ['#a60628', '#7a68a6', '#348abd'])):
    RSS = loadFile(filepath + spec['file'])
    o = RSS[246 - spec['fiber']]['object']
    e = RSS[246 - spec['fiber']]['error']
    o.smooth(15)
    e.smooth(15)
    if i == 1:
        ax[i].plot(o.wavelength / (1 + spec['z']), o.flux / 1.1, c=c)
        ax[i].fill_between(o.wavelength / (1 + spec['z']),
                           o.flux / 1.1 - e.flux,
                           e.flux + o.flux / 1.1,
                           color='0.8')
    else:
        ax[i].plot(o.wavelength / (1 + spec['z']), o.flux, c=c)
        ax[i].fill_between(o.wavelength / (1 + spec['z']),
                           o.flux - e.flux,
                           e.flux + o.flux,
                           color='0.8')

    ax[i].text(3250, 0.8, 'Q = %d' % i, ha='center', size=20)
    if i < 2:
        ax[i].text(4500, 0.2, 'z = %.4f' % spec['z'], ha='left', size=16)
    else:
        ax[i].text(4500, 0.2, 'z = ???', ha='left', size=16)

# plot some features
# k
ax[0].plot([3933.68, 3933.68], [0, 0.8], '--', c='k', lw=1)
ax[0].text(3933.68,
           0.95,
           'K',
           ha='center',
           va='top',
           size=11,
           rotation='vertical')
# H
ax[0].plot([3968.47, 3968.47], [0, 0.8], '--', c='k', lw=1)
ax[0].text(3968.47,
           0.95,
           'H',
           ha='center',
           va='top',
           size=11,
           rotation='vertical')
# G
ax[0].plot([4307.74, 4307.74], [0, 0.8], '--', c='k', lw=1)
ax[0].text(4307.74,
           0.95,
           'G',
           ha='center',
           va='top',
           size=11,
           rotation='vertical')

# [Oii]
ax[1].plot([3727.3, 3727.3], [0, 0.7], '--', c='k', lw=1)
ax[1].text(3727.3,
           0.95,
           '$[O_{II}]$',
           ha='center',
           va='top',
           size=11,
           rotation='vertical')

# tweak plot
for a in ax:
    a.set_xlim(3000, 5000)
    a.set_ylim(0, 1)
    a.set_yticks([0.2, 0.4, 0.6, 0.8])

ax[0].set_xticklabels([])
ax[1].set_xticklabels([])

ax[2].set_xlabel('Rest Wavelength (Angtrom)')
ax[1].set_ylabel('Relative Flux')

import pyfits as pyf
from glob import glob
import os
from astLib import astWCS
from numpy import array
import scipy

files = glob('bcs*.fits')

for f in files:
    oimg = pyf.open(f)
    obj = oimg[1].header['object'].split('_')
    print f, obj
    field, dither, num = obj[1].split()
    # Correct for a typo in the naming.
    if obj[0] == 'c205p08+46p7':
        obj[0] = 'c250p08+46p7'

    # load data and skysubtract
    data = oimg[1].data
    dataHDR = oimg[1].header
    odata_dim = data.shape
    wcs = astWCS.WCS(f, extensionName=1)
    owavelengthStartEnd = wcs.getImageMinMaxWCSCoords()[0:2]
    owavelengthStep = oimg[1].header['CDELT1']

    owavelengthRange = [owavelengthStartEnd[0] + i * owavelengthStep
                        for i in range(odata_dim[1])]

    if not os.path.isfile(obj[0]+'_'+field+'_'+num+'.fits'):
        # rewriting the whole file because that is easy to update
        oimg.writeto(obj[0]+'_'+field+'_'+num+'.fits')
        # update with sky subtraction
        #pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data, dataHDR, 1)
    else:
        # Here's the data we are going to add to
        img = pyf.open(obj[0]+'_'+field+'_'+num+'.fits')
        # Load the IFU data -- Row-stacked spectra
        odata = img[1].data
        oError = img[2].data
        odata_dim = odata.shape
        wcs = astWCS.WCS(obj[0]+'_'+field+'_'+num+'.fits', extensionName=1)
        owavelengthStartEnd = wcs.getImageMinMaxWCSCoords()[0:2]
        owavelengthStep = oimg[1].header['CDELT1']

        wavelengthRange = [owavelengthStartEnd[0] + i * owavelengthStep
                            for i in range(odata_dim[1])]


        data1 = array([scipy.interp(owavelengthRange, wavelengthRange, i) for i in
                odata])
        try:
            pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data1+data,
                    dataHDR, 1)
        except ValueError:
            print 'Different lengths'
            # Make sure all of the arrays are the same length
#            if data.shape[1] > data1.shape[1]:
                #sky.pop(-1*(data.shape[1]-data1.shape[1]))
#                data = delete(data, -1*(data.shape[1]-data1.shape[1]), 1)
#            elif data.shape[1] < data1.shape[1]:
#                data1 = delete(data1, -1*(data1.shape[1]-data.shape[1]), 1)
#            else:
#                print "I don't know what to do!"

            # UPDATE!!!
            pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data1+data,
                    dataHDR, 1)

        img.close()
    oimg.close()

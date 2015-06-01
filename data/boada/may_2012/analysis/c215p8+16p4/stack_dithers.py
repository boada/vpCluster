import pyfits as pyf
from glob import glob
import os
from astLib.astStats import clippedMedianStdev
from numpy import average, delete

import numpy as np

def skySubtract(data):
    sums = [sum(data[i,:]) for i in range(data.shape[0])]
    med = clippedMedianStdev(sums)
    med = med['clippedMedian']
    skyfibers = [i for i in range(data.shape[0]) if sum(data[i,:]) <= med]
    skydata = data.take(skyfibers, axis=0)
    skyflux = [average(skydata[:,i]) for i in range(skydata.shape[1])]

    return skyflux

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

    sky = skySubtract(data)

    if not os.path.isfile(obj[0]+'_'+field+'_'+num+'.fits'):
        # rewriting the whole file because that is easy to update
        oimg.writeto(obj[0]+'_'+field+'_'+num+'.fits')
        # update with sky subtraction
        pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data-sky, dataHDR, 1)
    else:
        # Here's the data we are going to add to
        img = pyf.open(obj[0]+'_'+field+'_'+num+'.fits')
        data1 = img[1].data
        dataHDR1 = img[1].header
        try:
            pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data1+data-sky,
                    dataHDR, 1)
        except ValueError:
            print 'Different lengths'
            # Make sure all of the arrays are the same length
            if data.shape[1] > data1.shape[1]:
                sky.pop(-1*(data.shape[1]-data1.shape[1]))
                data = delete(data, -1*(data.shape[1]-data1.shape[1]), 1)
            elif data.shape[1] < data1.shape[1]:
                data1 = delete(data1, -1*(data1.shape[1]-data.shape[1]), 1)
            else:
                print "I don't know what to do!"

            # UPDATE!!!
            pyf.update(obj[0]+'_'+field+'_'+num+'.fits', data1+data-sky,
                    dataHDR, 1)

        img.close()
    oimg.close()

import pyfits as pyf
import shutil
from glob import glob
import os

files = glob('*.fits')

for f in files:
    oimg = pyf.open(f)
    obj = oimg[1].header['object'].split('_')
    if not os.path.isdir(obj[0]):
        os.mkdir(obj[0])
    else:
        try:
            dummy = obj[1]
            print f, obj[0], obj[1]
        except IndexError:
            print f, obj[0]
#    shutil.copyfile(f, obj[0]+os.path.sep+f)

    oimg.close()

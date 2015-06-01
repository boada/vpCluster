import pyfits as pyf
import shutil
from glob import glob
import os

files = glob('unsorted/*.fits')

for f in files:
    oimg = pyf.open(f)
    obj = oimg[1].header['object'].split('_')
    if not os.path.isdir(obj[0]):
        os.mkdir(obj[0])
    else:
       print f, obj[0]
    shutil.copyfile(f, obj[0]+os.path.sep+f.split('/')[1])
    oimg.close()

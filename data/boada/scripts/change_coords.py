import os
from astLib import astCoords
import glob
import numpy as np

for files in glob.glob('*.txt'):
    fbits = files.split('_')
    if not os.path.isfile(fbits[0]+'_combined_tiles.txt'):
        tiles = glob.glob(fbits[0]+'*')
        with open(fbits[0]+'_combined_tiles.txt', 'wt') as newf:
            newf.writelines('# tile dither fiber fiber_ra fiber_dec\n')
            for f in tiles:
                tile = f.split('_')[1]
                dither = f.split('_')[2]
                print(fbits[0], tile, dither)
                data = np.genfromtxt(f, dtype='str')
                for fiber, ra, dec in zip(range(1,len(data[:,0])+1), data[:,1],
                    data[:,2]):
                    ra_deg = astCoords.hms2decimal(ra, ':')
                    dec_deg = astCoords.dms2decimal(dec, ':')
                    newf.writelines('%s %s %s %s %s\n' %
                        (tile, dither, fiber, ra_deg, dec_deg))




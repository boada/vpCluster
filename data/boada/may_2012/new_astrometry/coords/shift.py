import numpy as np
from astLib import astCoords
import sys

data = np.genfromtxt(sys.argv[1], dtype=str)

Ra_shift = 6.7995
Dec_shift = 1.4995

newf = open('shifted.txt', 'wt')
for fiber, ra, dec in zip(
        range(1, len(data[:, 0]) + 1), data[:, 1], data[:, 2]):
    ra_deg = astCoords.hms2decimal(ra, ':')
    dec_deg = astCoords.dms2decimal(dec, ':')
    new_coords = astCoords.shiftRADec(ra_deg, dec_deg, -1 * Ra_shift, -1 *
                                      Dec_shift)
    newf.writelines('%s %s %s\n' % (fiber, new_coords[0], new_coords[1]))

newf.close()

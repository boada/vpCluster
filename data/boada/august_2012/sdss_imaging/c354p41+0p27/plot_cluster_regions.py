import matplotlib
matplotlib.use('Agg')

import aplpy
import numpy as np
from astLib import astCoords

files = '/Users/steven/Projects/cluster/data/boada/august_2012/updated_astrometry_matched/coords/'

gc = aplpy.FITSFigure('fpC-002662-r4-0245.fit', north=True, figsize=(5, 5))
gc.show_grayscale(stretch='log')
gc.recenter(354.4155423, 0.2713716, radius=3 / 60.)
gc.show_circles(354.4155423, 0.2713716, radius=2.3 / 60.)

gc.set_auto_refresh(False)

fibers = np.loadtxt(files + 'c354p41+0p27_NE_D1_coords.txt', dtype='str')
for ra, dec in zip(fibers[:, 1], fibers[:, 2]):
    x = astCoords.hms2decimal(ra, ':')
    y = astCoords.dms2decimal(dec, ':')
    gc.show_circles(x, y, 2 / 3600., color='black')
fibers = np.loadtxt(files + 'c354p41+0p27_NW_D1_coords.txt', dtype='str')
for ra, dec in zip(fibers[:, 1], fibers[:, 2]):
    x = astCoords.hms2decimal(ra, ':')
    y = astCoords.dms2decimal(dec, ':')
    gc.show_circles(x, y, 2 / 3600., color='blue')
fibers = np.loadtxt(files + 'c354p41+0p27_SW_D1_coords.txt', dtype='str')
for ra, dec in zip(fibers[:, 1], fibers[:, 2]):
    x = astCoords.hms2decimal(ra, ':')
    y = astCoords.dms2decimal(dec, ':')
    gc.show_circles(x, y, 2 / 3600., color='green')
fibers = np.loadtxt(files + 'c354p41+0p27_SE_D1_coords.txt', dtype='str')
for ra, dec in zip(fibers[:, 1], fibers[:, 2]):
    x = astCoords.hms2decimal(ra, ':')
    y = astCoords.dms2decimal(dec, ':')
    gc.show_circles(x, y, 2 / 3600., color='purple')

gc.hide_tick_labels()
gc.hide_axis_labels()
gc.ticks.hide()
gc.set_theme('publication')

gc.refresh()

gc.save('pointing.eps')

import aplpy
import pylab as pyl

cluster = 'c203p83+41p00_r_mosaic.fits'

f = pyl.figure(1, figsize=(8,8))

gc = aplpy.FITSFigure(cluster, north=True, figure=f)
gc.show_grayscale(stretch='arcsinh')
gc.set_theme('publication')

gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm:ss')
#gc.set_tick_labels_size('small')



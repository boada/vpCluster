import aplpy
import pylab as pyl

cluster = 'c203p8+41p0_r_mosaic.fits'

f = pyl.figure(1, figsize=(8,8))

gc = aplpy.FITSFigure(cluster, north=True, figure=f)
gc.show_grayscale(stretch='arcsinh')
gc.set_theme('publication')

gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm:ss')
#gc.set_tick_labels_size('small')

data = pyl.genfromtxt(cluster.split('_')[0]+'_complete.csv', delimiter=',',
        names=True)

# filter out the specz's
x = pyl.isnan(data['Specz'])

# draw the specz's
gc.show_markers(data['ra'][~x], data['dec'][~x], edgecolor='#e24a33',
        facecolor='none', marker='D', s=200)



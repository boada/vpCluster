import aplpy
import pylab as pyl
import sys

#cluster = 'c260p61+32p13_r_mosaic.fits'
def main(cluster):
    cluster = cluster + '_r_mosaic.fits'
    f = pyl.figure(1, figsize=(8,8))

    gc = aplpy.FITSFigure(cluster, north=True, figure=f)
    gc.show_grayscale(stretch='arcsinh')
    gc.set_theme('publication')

    gc.set_tick_labels_format(xformat='hh:mm:ss', yformat='dd:mm:ss')
    #gc.set_tick_labels_size('small')

    data = pyl.genfromtxt('./../analysis_all/' +\
            cluster.split('_')[0]+'_redshifts.csv', delimiter=',', names=True,
            dtype=None)

    try:
        # filter out the specz's
        x = pyl.isnan(data['Specz'])
        # draw the specz's
        gc.show_markers(data['ra'][~x], data['dec'][~x], edgecolor='#e24a33',
                facecolor='none', marker='D', s=200)

    except ValueError:
        print 'no Speczs found'

    # draw observed but not redshifted
    x = (~pyl.isnan(data['fiber'])) & (pyl.isnan(data['redshift']))
    gc.show_markers(data['ra'][x], data['dec'][x], edgecolor='#a60628',
            facecolor='none', marker='s', s=150)

    # draw redshifted
    x = (~pyl.isnan(data['fiber'])) & (~pyl.isnan(data['redshift']))
    gc.show_markers(data['ra'][x], data['dec'][x], edgecolor='#188487',
            facecolor='none', marker='o', s=150)
if __name__ == "__main__":
    main(sys.argv[1])

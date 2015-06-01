import pandas as pd
cluster = './c328p33+0p19'
photoz = pd.read_csv(cluster+'_photoz.csv')
specz = pd.read_csv(cluster+'_specz.csv')
photo = pd.read_csv(cluster+'_photo.csv')
c = pd.merge(photo, photoz, left_on='objID', right_on='objID', how='outer')
try:
    c = pd.merge(c, specz, left_on='objID', right_on='objID', how='outer')
    sc.rename(columns={'ra_x':'ra', 'dec_x':'dec'}, inplace=True)
    del sc['ra']
    del sc['dec']
    del sc['dec_y']
    del sc['ra_y']
    sc.to_csv(cluster+'_combined.csv')
    matched = pd.read_csv(cluster+'_matched.csv')
    scm = pd.merge(matched, sc, left_on='objID', right_on='objID', how='outer')
except KeyError:
    print('No Speczs')
    c.to_csv(cluster+'_combined.csv')
    matched = pd.read_csv(cluster+'_matched.csv')
    scm = pd.merge(matched, c, left_on='objID', right_on='objID', how='outer')

scm.to_csv(cluster+'_complete.csv')

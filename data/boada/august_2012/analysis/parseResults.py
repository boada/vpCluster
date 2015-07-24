import glob
from numpy import genfromtxt, asarray, savetxt

files = glob.glob('*_results/*.results')

r = []
for f in files:
    print f
    cluster, field, dither = f.split('/')[1].split('_')
    data = genfromtxt(f, delimiter='\t', names=True, dtype=None)
    try:
        for fiber, z, Q in zip(data['Fiber'], data['Redshift'],
                data['Quality']):
            if Q == 0:
                r.append((cluster, field, dither.split('.')[0], fiber, z))
    except TypeError:
        fiber = int(data['Fiber'])
        z = float(data['Redshift'])
        Q = int(data['Quality'])
        if Q == 0:
            r.append((cluster, field, dither.split('.')[0], fiber, z))

print len(r), 'objects read'
#print r

r = asarray(r)
#data = zeros((r.shape[0],), dtype=[('cluster', 'a13'), ('field', 'a2'), ('dither',
#    '>i4'), ('fiber', '>i4'), ('redshift', '>f4')])
#data['cluster'] = r[:,0]
#data['field'] = r[:,1]
#data['dither'] = r[:,2]
#data['fiber'] = r[:,3]
#data['redshift'] = r[:,4]

savetxt('results_august12', r, delimiter=' ', fmt='%s')

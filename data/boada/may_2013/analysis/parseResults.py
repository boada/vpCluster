import pylab as pyl
import glob
from numpy import genfromtxt

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
                r.append((field, dither.split('.')[0], fiber, z))
    except TypeError:
        fiber = int(data['Fiber'])
        z = float(data['Redshift'])
        Q = int(data['Quality'])
        if Q == 0:
            r.append((field, dither, fiber, z))

print len(r), 'objects read'
print r

print len(r)
#pyl.hist(r)
#pyl.show()

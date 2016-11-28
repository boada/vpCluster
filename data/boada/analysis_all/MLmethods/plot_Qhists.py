import h5py as hdf
import pylab as pyl

for i in range(20):
    data_path = '../data/buzzard_v1.0/allbands/truth'
    with hdf.File(data_path+'/truth'+str(i).zfill(2)+'_Oii.hdf5', 'r') as f:

        mags = f['truth%s_Oii' % (str(i).zfill(2))]
        qs = f['Q']
        try:
            rMag = pyl.append(rMag, mags['OMAG'][:,2]) # r band
            Qs = pyl.append(Qs, qs.value)
        except NameError:
            rMag = mags['OMAG'][:,2]
            Qs = qs.value

q0 = pyl.where(Qs == 0)[0]
r0 = rMag[q0]
q1 = pyl.where(Qs == 1)[0]
r1 = rMag[q1]
q2 = pyl.where(Qs == 2)[0]
r2 = rMag[q2]

# make a figure
f = pyl.figure(figsize=(5,5*(pyl.sqrt(5.)-1.0)/2.0))
ax = f.add_subplot(111)

bins = pyl.linspace(14,22,15)
ax.hist(r2, weights=pyl.zeros_like(r2)+1./r2.size, histtype='step', bins=bins,
        lw=2, label='Q=2')
ax.hist(r1, weights=pyl.zeros_like(r1)+1./r1.size, histtype='step', bins=bins,
        lw=2, label='Q=1')
ax.hist(r0, weights=pyl.zeros_like(r0)+1./r0.size, histtype='step', bins=bins,
        lw=2, label='Q=0')

ax.legend(loc='upper right')
ax.invert_xaxis()

ax.set_ylim(0,0.5)
ax.set_xlabel('$m_r$ (mag)')
ax.set_ylabel('Fraction of Total')
pyl.show()


# final fractions
tot = pyl.where(Qs != -1)[0]
print q0.size/float(tot.size)
print q1.size/float(tot.size)
print q2.size/float(tot.size)


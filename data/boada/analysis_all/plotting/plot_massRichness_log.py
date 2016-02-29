from __future__ import division # force float division
import pylab as pyl

def mkError(x, x_err):

    return 0.434 * x_err/x

richnessData = pyl.genfromtxt('./boada_rich.txt', names=True, dtype=None)
massData = pyl.genfromtxt('./cluster_props', names=True, dtype=None)

richnessData.sort()
massData.sort()

yerr = [mkError(massData['mass'], massData['mass']-massData['mass_lower']),
        mkError(massData['mass'], abs(massData['mass'] -\
            massData['mass_upper']))]


#yerr = [(massData['mass']-massData['mass_lower'])*10, abs(massData['mass'] -
#    massData['mass_upper'])*10]

ax = pyl.subplot()

ax.errorbar(pyl.log10(richnessData['lambda']), pyl.log10(massData['mass']),
        xerr=mkError(richnessData['lambda'], richnessData['lambda_err']),
        yerr=yerr, fmt='o', label='This Work')

ax.set_xlabel('Log Richness')
ax.set_ylabel(r'Log Cluster Mass $[M_{\odot}]$')

# add the Rozo2010 points

rozo = pyl.genfromtxt('./rozo2010_points', names=True, dtype=None)
ax.errorbar(pyl.log10(rozo['lambda']), pyl.log10(rozo['mass']*1e15),
        xerr=mkError(rozo['lambda'], rozo['lambda_err']), fmt='o',
        label='Rozo+2010')

pyl.legend(loc='lower right')

#ax.set_ylim(1e-2, 50)
#pyl.semilogy()
pyl.show()


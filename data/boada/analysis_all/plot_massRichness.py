import pylab as pyl

richnessData = pyl.genfromtxt('./boada_rich.txt', names=True, dtype=None)
massData = pyl.genfromtxt('./cluster_props', names=True, dtype=None)

richnessData.sort()
massData.sort()

yerr = [(massData['mass']-massData['mass_lower'])*10, abs(massData['mass'] -
    massData['mass_upper'])*10]

ax = pyl.subplot()

ax.errorbar(richnessData['lambda'], massData['mass']*10,
        xerr=richnessData['lambda_err'], yerr=yerr, fmt='o', label='This Work')

ax.set_xlabel('Richness')
ax.set_ylabel(r'Cluster Mass $[10^{14}M_{\odot}]$')

# add the Rozo2010 points

rozo = pyl.genfromtxt('./rozo2010_points', names=True, dtype=None)
ax.errorbar(rozo['lambda'], rozo['mass']*10,
        xerr=rozo['lambda_err'], fmt='o', label='Rozo+2010')

pyl.legend(loc='lower right')

ax.set_ylim(1e-2, 50)
pyl.semilogy()
pyl.show()


import pandas as pd
from glob import glob
from calc_mass import calc_mass
from mod_dominance import mod_dom
import pylab as pyl
import matplotlib.mlab as mlab

files = glob('members/*.csv')

i = 1
for f in files:
    data = pd.read_csv(f)
    data, mass, vd = calc_mass(data)
    mod, dom = mod_dom(f)

    pyl.figure(i)
    result = pyl.hist(data['LOSV'])
    x = pyl.linspace(data['LOSV'].min(), data['LOSV'].max(), 100)
    dx = result[1][1] - result[1][0]
    scale = len(data['LOSV'])*dx
    pyl.plot(x, mlab.normpdf(x, 0, vd) * scale, lw=2, c='k')
    pyl.xlabel('Line of Sight Velocity (km/s)')
    pyl.ylabel('Count')
    pyl.title('Mod = %s Dom = %s' %(mod, dom))
    pyl.savefig(str(i)+'.png')
    i+=1

pyl.show()

import sys
import numpy as np


def main(filein):

    data = np.genfromtxt(filein, names=True, delimiter=',')

    f1 = open(filein.rstrip('csv') + 'reg', 'wt')
    f1.writelines('# Region file format: DS9 version 4.1\n')
    f1.writelines('# Filename:' + filein.rstrip('csv') + 'reg\n')

    for ra, dec, in zip(data['ra'], data['dec']):
        f1.writelines('fk5;circle(' + str(ra) + ',')
        f1.writelines(str(dec) + ',2")\n')

    f1.close()


if __name__ == "__main__":
    main(sys.argv[1])

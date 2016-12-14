import pylab as pyl


def compare(features, targets, columns):
    q0mask = pyl.where(targets == 0)[0]
    q1mask = pyl.where(targets == 1)[0]
    q2mask = pyl.where(targets == 2)[0]

    pyl.scatter(features[:, columns[0]][q0mask],
                features[:, columns[1]][q0mask], c='r', label='0')
    pyl.scatter(features[:, columns[0]][q1mask],
                features[:, columns[1]][q1mask], c='g', label='1')
    pyl.scatter(features[:, columns[0]][q2mask],
                features[:, columns[1]][q2mask], c='b', label='2')

    pyl.legend(loc='best')
    pyl.show()

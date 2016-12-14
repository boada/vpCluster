#!/home/users/boada/resources/python-2.7.5/bin/python

import urllib
import sys
import os
import string
from astLib import astCoords
import numpy as np
import glob
import time


def filtercomment(sql):
    "Get rid of comments starting with --"
    fsql = ''
    for line in sql.split('\n'):
        fsql += line.split('--')[0] + ' ' + os.linesep
    return fsql


def query(sql):
    url = 'http://cas.sdss.org/public/en/tools/search/x_sql.asp'
    fmt = 'csv'
    "Run query and return file object"
    fsql = filtercomment(sql)
    params = urllib.urlencode({'cmd': fsql, 'format': fmt})
    return urllib.urlopen(url + '?%s' % params)


def work(ra, dec, outfile):
    select = '''SELECT TOP 100
    G.objID,
    G.ra,
    G.dec,
    G.u,
    G.Err_u,
    G.g,
    G.Err_g,
    G.r,
    G.Err_r,
    G.i,
    G.Err_i,
    G.z,
    G.Err_z,
    Pz.z,
    Pz.zErr

    '''
    FROM = '''FROM
    Galaxy as G

    '''
    join = '''JOIN
    dbo.fGetNearbyObjEq(''' + str(ra) + ',' + str(dec) + ''',0.033) AS GN
    ON
    G.objID = GN.objID

    JOIN
    Photoz as Pz
    ON
    G.objID = Pz.objID

    '''
    where = '''Where G.g < 22'''

    sql = select + FROM + join + where

    result = query(sql)

    line = result.readline()
    if line.startswith("ERROR"):
        ofp = sys.stderr
    else:
        ofp = open(outfile, 'wt')
        while line:
            ofp.write(string.rstrip(line) + os.linesep)
            line = result.readline()
        ofp.close()

# get file list
filelist = glob.glob('*D1*') + glob.glob('*D2*') + glob.glob('*D3*')

cwd = os.getcwd()

for f in filelist:
    data = np.loadtxt(f, dtype='str')
    fBits = f.split('_')
    print fBits[0]
    if not os.path.isdir(fBits[0]):
        os.mkdir(fBits[0])
    os.chdir(fBits[0])
    if not os.path.isdir(fBits[1]):
        os.mkdir(fBits[1])
    os.chdir(fBits[1])
    if not os.path.isdir(fBits[2]):
        os.mkdir(fBits[2])
    os.chdir(fBits[2])

    for fiber, ra, dec in zip(data[:, 0], data[:, 1], data[:, 2]):
        if not os.path.isfile(fiber.zfill(3) + '.txt'):
            ra = astCoords.hms2decimal(ra, ':')
            dec = astCoords.dms2decimal(dec, ':')
            work(ra, dec, fiber.zfill(3) + '.txt')
            #print fiber, ra, dec

            time.sleep(2)
        else:
            pass

    os.chdir(cwd)

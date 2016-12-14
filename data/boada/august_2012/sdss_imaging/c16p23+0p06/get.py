import glob
import subprocess

filelist = glob.glob('sdss*')

for f in filelist:
    ofp = open(f, 'rt')
    line = ofp.readline()
    while line:
        p = subprocess.Popen('curl -O ' + line,
                             shell=True,
                             stdout=subprocess.PIPE)
        p.wait()
        line = ofp.readline()

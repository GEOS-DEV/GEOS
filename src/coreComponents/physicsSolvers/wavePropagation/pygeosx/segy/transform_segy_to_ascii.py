import sys
from obspy.io.segy.segy import _read_segy
from segyio import open
import numpy as np

filename=str(sys.argv[1])
#segy = _read_segy(filename)
segy = open(filename)
# turn ObsPy Stream in a matrix of traces
# first dimension time, second dimension traces
ntraces = len(segy.trace)
nsamples = len(segy.trace[0].data)
mtraces = np.zeros((nsamples, ntraces))
i = 0
for tr in segy.trace:
    mtraces[:, i] = tr.data[:]
    i += 1

np.savetxt("output.txt",mtraces,fmt="%s")



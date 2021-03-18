import sys
from obspy.io.segy.segy import _read_segy
import numpy as np
import matplotlib.pyplot as plt


if (len(sys.argv)==4):
    filename=str(sys.argv[1])
    filename2=str(sys.argv[2])
    nb_trace=int(sys.argv[3])+1
    segy = _read_segy(filename)
    segy2 = _read_segy(filename2)

    # turn ObsPy Stream in a matrix of traces
    # first dimension time, second dimension traces
    ntraces = len(segy.traces)
    nsamples = len(segy.traces[0].data)
    mtraces = np.zeros((nsamples, ntraces))
    i = 0
    for tr in segy.traces:
        mtraces[:, i] = tr.data[:]
        i += 1

    # turn ObsPy Stream in a matrix of traces
    # first dimension time, second dimension traces
    ntraces2 = len(segy2.traces)
    nsamples2 = len(segy2.traces[0].data)
    mtraces2 = np.zeros((nsamples2, ntraces2))
    i = 0
    for tr in segy2.traces:
        mtraces2[:, i] = tr.data[:]
        i += 1

        
    graph=plt.plot(range(nsamples),mtraces[:,nb_trace]-mtraces2[:,nb_trace],label='graph',color='blue')
    plt.xlabel("time samples")
    plt.ylabel("pressure")
    plt.legend(["Diff Trace n°"+str(nb_trace-1)])
    
    plt.show()
else :
    print('Error, argument(s) invalid')

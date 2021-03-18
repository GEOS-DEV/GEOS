import sys
from obspy.io.segy.segy import _read_segy
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv)==3):
    filename=str(sys.argv[1])
    filename2=str(sys.argv[2])
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


    fig, ax = plt.subplots()
    #ax.axis('equal')

    error = 0.0
    error = sum(sum((mtraces-mtraces2)**2))
    print('ERROR :',error)

    plt.imshow(mtraces-mtraces2, aspect='auto', interpolation='none')
    fig.tight_layout()
    plt.set_cmap('seismic')
    plt.colorbar()
    plt.show()


elif (len(sys.argv)==4):
    filename=str(sys.argv[1])
    filename2=str(sys.argv[2])
    min_val=-abs(float((sys.argv[3])))
    max_val=abs(float((sys.argv[3])))
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

    error = 0.0
    error = sum(sum((mtraces-mtraces2)**2))
    print('ERROR :',error)


    fig, ax = plt.subplots()
    #ax.axis('equal')

    plt.imshow(mtraces-mtraces2, aspect='auto', interpolation='none')
    fig.tight_layout()
    plt.clim(vmin=min_val,vmax=max_val)
    plt.set_cmap('seismic')
    plt.colorbar()
    plt.show()


else :
    print('Error, argument(s) invalid')

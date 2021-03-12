import sys
import numpy as np
import matplotlib.pyplot as plt

if (len(sys.argv)==2):
    filename=str(sys.argv[1])

    #ntraces = len(segy.traces)
    #nsamples = len(segy.traces[0].data)
    #mtraces = np.zeros((nsamples, ntraces))
    #i = 0
    #for tr in segy.traces:
    #    mtraces[:, i] = tr.data[:]
    #    i += 1

    ####Read ASCII file####
    mtraces=np.loadtxt(filename)

    fig, ax = plt.subplots()
    #ax.axis('equal')

    plt.imshow(mtraces, aspect='auto', interpolation='none')
    fig.tight_layout()
#    plt.set_cmap('seismic')
    plt.set_cmap('Greys')
    plt.colorbar()
    plt.show()


elif (len(sys.argv)==3):
    filename=str(sys.argv[1])
    min_val=-abs(float((sys.argv[2])))
    max_val=abs(float((sys.argv[2])))
    #segy = _read_segy(filename)

    ## turn ObsPy Stream in a matrix of traces
    ## first dimension time, second dimension traces
    #ntraces = len(segy.traces)
    #nsamples = len(segy.traces[0].data)
    #mtraces = np.zeros((nsamples, ntraces))
    #i = 0
    #for tr in segy.traces:
    #    mtraces[:, i] = tr.data[:]
    #    i += 1

    ####Read ASCII file####
    mtraces=np.loadtxt(filename)

    fig, ax = plt.subplots()
    #ax.axis('equal')

    plt.imshow(mtraces, aspect='auto', interpolation='none')
    fig.tight_layout()
    plt.clim(vmin=min_val,vmax=max_val)
#    plt.set_cmap('seismic')
    plt.set_cmap('Greys')
    plt.colorbar()
    plt.show()


elif (len(sys.argv)==4):
    filename=str(sys.argv[1])
    min_val=float((sys.argv[2]))
    max_val=float((sys.argv[3]))
    #segy = _read_segy(filename)

    ## turn ObsPy Stream in a matrix of traces
    ## first dimension time, second dimension traces
    #ntraces = len(segy.traces)
    #nsamples = len(segy.traces[0].data)
    #mtraces = np.zeros((nsamples, ntraces))
    #i = 0
    #for tr in segy.traces:
    #    mtraces[:, i] = tr.data[:]
    #    i += 1

    ####Read ASCII file####
    mtraces=np.loadtxt(filename)

    fig, ax = plt.subplots()
    #ax.axis('equal')

    plt.imshow(mtraces, aspect='auto', interpolation='none')
    fig.tight_layout()
    plt.clim(vmin=min_val,vmax=max_val)
    plt.set_cmap('seismic')
    plt.colorbar()
    plt.show()

else :
    print('Error, argument(s) invalid')

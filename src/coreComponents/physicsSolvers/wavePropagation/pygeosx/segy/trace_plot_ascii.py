import sys
import numpy as np
import matplotlib.pyplot as plt


if (len(sys.argv)==3):
    filename=str(sys.argv[1])
    nb_trace=int(sys.argv[2])+1
    #segy = _read_segy(filename)

    # turn ObsPy Stream in a matrix of traces
    # first dimension time, second dimension traces
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
    plt.gca().invert_yaxis()
    graph=plt.plot(mtraces[:,nb_trace],range(nsamples),label='graph',color='blue')
    plt.xlabel("time samples")
    plt.ylabel("pressure")
    plt.legend(["Trace n°"+str(nb_trace-1)])

    plt.show()


elif (len(sys.argv)==4):

    filename=str(sys.argv[1])
    nb_trace=int(sys.argv[2])+1
    min_val=-abs(float((sys.argv[3])))
    max_val=abs(float((sys.argv[3])))
    #segy = _read_segy(filename)

    # turn ObsPy Stream in a matrix of traces
    # first dimension time, second dimension traces
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
    plt.gca().invert_yaxis()
    plt.axis([min_val, max_val, nsamples, 1])
    graph=plt.plot(mtraces[:,nb_trace],range(nsamples),label='graph',color='blue')
    plt.xlabel("time samples")
    plt.ylabel("pressure")
    plt.legend(["Trace n°"+str(nb_trace-1)])



    fichier = open("data.txt", "w")
    for i in range(nsamples):
        fichier.write(str(i)+" "+str(mtraces[i,nb_trace])+"\n")

    fichier.close()
    plt.show()


else :
    print('Error, argument(s) invalid')

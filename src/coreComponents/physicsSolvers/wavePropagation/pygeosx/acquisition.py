'''
Created on 9/02/2021

@author: macpro
'''
import numpy as np
import copy
import segyio
import random
import os, glob

from shot import *
from receiver import *
from source import *
from segyManager import export_for_acquisition

from mpmain import segyPath


def segy_acquisition(folder_path, wavelet, dt):

    shot_list = []

    for filename in glob.glob(os.path.join(folder_path, '*.sgy')):
        receiver_list = []

        with segyio.open(filename, 'r', ignore_geometry=True) as f:
            scalarXY = float(f.header[0][71])
            scalarZ = float(f.header[0][69])

            sourceX = f.header[0][73]*abs(scalarXY)**np.sign(scalarXY)
            sourceY = f.header[0][77]*abs(scalarXY)**np.sign(scalarXY)
            sourceZ = f.header[0][49]*abs(scalarXY)**np.sign(scalarZ)

            source = Source([sourceX, sourceY, sourceZ], wavelet, dt)
            for i in range(len(f.trace)):

                receiverX = f.header[i][81]*abs(scalarXY)**np.sign(scalarXY)
                receiverY = f.header[i][85]*abs(scalarXY)**np.sign(scalarXY)
                receiverZ = f.header[i][41]*abs(scalarXY)**np.sign(scalarZ)

                receiver_list.append(Receiver([receiverX, receiverY, receiverZ]))

            receivers = ReceiverSet(receiver_list)

            shot_list.append(Shot(source, receivers))

    return shot_list



def moving_acquisition(box,
                       wavelet,
                       dt,
                       nbsourcesx = 1,
                       nbsourcesy = 1,
                       nbreceiversx = 1,
                       nbreceiversy = 1,
                       lenRx = 0,
                       lenRy = 0,
                       export = 0):

    """Marine seismic acquisition reprensenting a boat pulling a source
       and a set of receivers in a squared domain

    Parameters
    ----------
    box :
        Numpy array containing min/max boundary coordinates of the domain

    wavelet :
        Source function (Ricker)

    nbsourcesx :
        Number of sources along x axis

    nbsourcesy :
        Number of sources along y axis

    nbreceiversx :
        Number of sources along x axis

    nbreceiversy :
        Number of sources along y axis

    lenRx :
        Lenght of the area for receivers along x axis

    lenRy :
        Lenght of the area for receivers along y axis

    Return
    ------
    shots :
        A list of Shot objects

    Notes
    -----
    - "eps" variable is defined here to avoid source and receivers
      to be on a node coord or on an element edge
    - This is a temporary acquisition for simulation tests
    - Seismic acquisition will be read from segy files
    """

    eps  = 49
    xmin = box[0][0]
    xmax = box[0][1]

    ymin = box[1][0]
    ymax = box[1][1]

    zmin = box[2][0]
    zmax = box[2][1] - eps

    movex=(xmax-xmin)/(nbsourcesx+1) + 1
    movey=(ymax-ymin)/(nbsourcesy+1) + 1

    srcpos=[xmin,ymin,zmax]

    xposleft = np.linspace(xmin + 0.3*movex, xmin + 0.3*movex + lenRx, nbreceiversx)
    xposright = np.linspace(xmax - 0.3*movex, xmax - 0.3*movex - lenRx, nbreceiversx)

    ypos = np.linspace(ymin - lenRy/2, ymin + lenRy/2, nbreceiversy)

    receiversbaseleft = ReceiverSet([Receiver([x, y, zmax]) for x in xposleft for y in ypos])
    receiversbaseright = ReceiverSet([Receiver([x, y, zmax]) for x in xposright for y in ypos])

    shots=[]

    for i in range(nbsourcesy):
        if i%2==0:
            srcpos[0] = xmin
            receivers = copy.deepcopy(receiversbaseleft)

        else:
            srcpos[0] = xmax
            receivers = copy.deepcopy(receiversbaseright)

        srcpos[1] += movey
        receivers.linearSetMove(1, (i+1)*movey)
        for j in range(nbsourcesx):

            if i%2==0:
                srcpos[0] += movex
                receivers.linearSetMove(0, movex)

            else:
                srcpos[0] -= movex
                receivers.linearSetMove(0, -movex)

            rcvset = receivers.getInsideDomain(box)

            # Define source location and type
            source = Source(srcpos, wavelet, dt)

            shot = Shot(source, copy.deepcopy(rcvset))
            shots.append(shot)

    if export:
        acq_name = "/moving_Sx=" + str(nbsourcesx) +"_Sy=" + str(nbsourcesy) + "_Rx=" + str(nbreceiversx) + "_Ry=" + str(nbreceiversy) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots



def random_acquisition(box,
                       wavelet,
                       dt,
                       nbsources,
                       nbreceiversx,
                       nbreceiversy,
                       export = 0):
    """Random seismic acquisition, the positions of sources are set randomly,
       the receivers are set as a grid over the domain

    Parameters
    ----------
    box :
        Numpy array containing min/max boundary coordinates of the domain

    wavelet :
        Source function (Ricker)

    nbsources :
        Number of sources

    nbreceiversx :
        Number of sources along x axis

    nbreceiversy :
        Number of sources along y axis

    Return
    ------
    shots :
        A list of Shot objects
    """

    xmin = box[0][0]
    xmax = box[0][1]

    ymin = box[1][0]
    ymax = box[1][1]

    zmin = box[2][0]
    zmax = box[2][1]


    xpos = np.linspace(xmin, xmax, nbreceiversx)
    ypos = np.linspace(ymin, ymax, nbreceiversy)

    receivers = ReceiverSet([Receiver([x, y, zmax]) for x in xpos for y in ypos])

    shots=list()

    for i in range(nbsources):

        srcpos = [random.random()*(xmin-0.2*(xmax-xmin)),(xmax+0.2*(xmax-xmin))*random.random(), zmax]

        # Define source location and type
        source = Source(srcpos, wavelet, dt)

        shot = Shot(source, receivers)
        shots.append(shot)

    if export:
        acq_name = "/random_nbS=" + str(nbsources) + "_Rx=" + str(nbreceiversx) + "_Ry=" + str(nbreceiversy) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots



'''Create equiscaped acquisition (copy from Pysit)'''
def equispaced_acquisition(box,
                           wavelet,
                           dt,
                           sourcex,
                           sourcey,
                           depthSource,
                           receiversx,
                           receiversy,
                           depthReceivers,
                           export = 0
                           ):

    xr = np.linspace(receiversx[0], receiversx[1], receiversx[2])
    yr = np.linspace(receiversy[0], receiversy[1], receiversy[2])

    xs = np.linspace(sourcex[0], sourcex[1], sourcex[2])
    ys = np.linspace(sourcey[0], sourcey[1], sourcey[2])

    receivers = ReceiverSet([Receiver([x, y, depthReceivers]) for x in xr for y in yr])

   # receivers = receivers.getInsideDomain(box)

    shots = []

    for i in range(len(ys)):
        for j in range(len(xs)):

            srcpos = [xs[j], ys[i], depthSource]
            source = Source(srcpos, wavelet, dt)

            shot = Shot(source, receivers)
            shots.append(shot)

    if export:
        acq_name = "/equispaced_Sx=" + str(len(xs)) +"_Sy=" + str(len(ys)) + "_Rx=" + str(len(xr)) + "_Ry=" + str(len(yr)) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots

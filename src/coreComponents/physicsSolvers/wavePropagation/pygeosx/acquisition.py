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



def moving_acquisition(boundaries,
                       wavelet,
                       dt,
                       number_of_sources_x = 1,
                       number_of_sources_y = 1,
                       number_of_receivers_x = 1,
                       number_of_receivers_y = 1,
                       source_depth = 1,
                       receivers_depth = 1,
                       zone_receiver_x = 0,
                       zone_receiver_y = 0,
                       export = 0,
                       ):


    """Marine seismic acquisition reprensenting a boat pulling a source
       and a set of receivers in a squared domain

    Parameters
    ----------
    box :
        List containing min/max boundary coordinates of the domain [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

    wavelet :
        Source function (Ricker)

    number_of_sources_x :
        Number of sources along x axis (this is the number of configuration, 1 position for each configuration)

    number_of_sources_y :
        Number of sources along y axis

    number_of_receivers_x :
        Number of sources along x axis

    number_of_receivers_y :
        Number of sources along y axis

    zone_receiver_x :
        Lenght of the area for receivers along x axis (distance from the receivers coordinate)

    zone_receiver_y :
        Lenght of the area for receivers along y axis

    Return
    ------
    shots :
        A list of Shot objects
    """


    xmin = boundaries[0][0]
    xmax = boundaries[0][1]

    ymin = boundaries[1][0]
    ymax = boundaries[1][1]

    zmin = boundaries[2][0]
    zmax = boundaries[2][1]

    movex=(xmax-xmin)/(number_of_sources_x + 1) + 1
    movey=(ymax-ymin)/(number_of_sources_y + 1) + 1

    srcpos=[xmin, ymin, source_depth]

    xposleft = np.linspace(xmin + 0.3*movex, xmin + 0.3*movex + zone_receiver_x, number_of_receivers_x)
    xposright = np.linspace(xmax - 0.3*movex, xmax - 0.3*movex - zone_receiver_x, number_of_receivers_x)

    ypos = np.linspace(ymin - zone_receiver_y/2, ymin + zone_receiver_y/2, number_of_receivers_y)

    receiversbaseleft = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xposleft for y in ypos])
    receiversbaseright = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xposright for y in ypos])

    shots=[]

    for i in range(number_of_sources_y):
        if i%2==0:
            srcpos[0] = xmin
            receivers = copy.deepcopy(receiversbaseleft)

        else:
            srcpos[0] = xmax
            receivers = copy.deepcopy(receiversbaseright)

        srcpos[1] += movey
        receivers.linearSetMove(1, (i+1)*movey)
        for j in range(number_of_sources_x):

            if i%2==0:
                srcpos[0] += movex
                receivers.linearSetMove(0, movex)

            else:
                srcpos[0] -= movex
                receivers.linearSetMove(0, -movex)

            rcvset = receivers.getInsideDomain(boundaries)

            # Define source location and type
            source = Source(srcpos, wavelet, dt)

            shot = Shot(source, copy.deepcopy(rcvset))
            shots.append(shot)

    if export:
        acq_name = "/moving_SxSy=" + str(number_of_sources_x) + "x" + str(number_of_sources_y) + "_RxRy=" + str(number_of_receivers_x) + "x" + str(number_of_receivers_y) + "_zoneR=" + str(2*zone_receiver_x) + "x" + str(2*zone_receiver_y) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots



def random_acquisition(boundaries,
                       wavelet,
                       dt,
                       number_of_sources = 1,
                       source_depth = 1,
                       number_of_receivers_= 1,
                       number_of_receivers_y = 1,
                       receivers_depth = 1,
                       export = 0):
    """Random seismic acquisition, the positions of sources are set randomly,
       the receivers are set as a grid over the domain

    Parameters
    ----------
    box :
        List containing min/max boundary coordinates of the domain [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

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

    xmin = boundaries[0][0]
    xmax = boundaries[0][1]

    ymin = boundaries[1][0]
    ymax = boundaries[1][1]

    zmin = boundaries[2][0]
    zmax = boundaries[2][1]


    xpos = np.linspace(xmin, xmax, number_of_receivers_x)
    ypos = np.linspace(ymin, ymax, number_of_receivers_y)

    receivers = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xpos for y in ypos])

    shots = list()

    for i in range(number_of_sources):

        srcpos = [random.random()*(xmin-0.2*(xmax-xmin)),(xmax+0.2*(xmax-xmin))*random.random(), source_depth]

        # Define source location and type
        source = Source(srcpos, wavelet, dt)

        shot = Shot(source, receivers)
        shots.append(shot)

    if export:
        acq_name = "/random_nbS=" + str(number_of_sources) + "_RxRy=" + str(number_of_receivers_x) + "x" + str(number_of_receivers_y) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots



'''Create equiscaped acquisition (copy from Pysit)'''
def equispaced_acquisition(boundaries,
                           wavelet,
                           dt,
                           source_boundaries_pos_x,
                           source_boundaries_pos_y,
                           receivers_boundaries_pos_x,
                           receivers_boundaries_pos_y,
                           number_of_sources_x,
                           number_of_sources_y,
                           number_of_receivers_x,
                           number_of_receivers_y,
                           source_depth = 1,
                           receivers_depth = 1,
                           export = 0
                           ):

    xr = np.linspace(receivers_boundaries_pos_x[0], receivers_boundaries_pos_x[1], number_of_receivers_x)
    yr = np.linspace(receivers_boundaries_pos_y[0], receivers_boundaries_pos_y[1], number_of_receivers_y)

    xs = np.linspace(source_boundaries_pos_x[0], source_boundaries_pos_x[1], number_of_sources_x)
    ys = np.linspace(source_boundaries_pos_y[0], source_boundaries_pos_y[1], number_of_sources_y)

    receivers = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xr for y in yr])

   # receivers = receivers.getInsideDomain(box)

    shots = []

    for i in range(len(ys)):
        for j in range(len(xs)):

            srcpos = [xs[j], ys[i], source_depth]
            source = Source(srcpos, wavelet, dt)

            shot = Shot(source, receivers)
            shots.append(shot)

    if export:
        acq_name = "/equispaced_SxSy=" + str(len(xs)) +"x" + str(len(ys)) + "_RxRy=" + str(len(xr)) + "x" + str(len(yr)) + "/"
        export_for_acquisition(shots, segyPath, acq_name)

    return shots


def cross_acquisition(box,
                      wavelet,
                      dt,
                      source_boundaries_pos_x,
                      source_boundaries_pos_y,
                      receivers_boundaries_pos_x,
                      receivers_boundaries_pos_y,
                      number_of_sources_x = 1,
                      number_of_sources_y = 1,
                      source_depth = 1,
                      number_of_receivers_x = 1,
                      number_of_receivers_y = 1,
                      receivers_depth = 1,
                      export = 0
                      ):


    xr1 = np.linspace(receivers_boundaries_pos_x[0][0], receivers_boundaries_pos_x[0][1], number_of_receivers_x)
    yr1 = np.linspace(receivers_boundaries_pos_y[0][0], receivers_boundaries_pos_y[0][1], number_of_receivers_y)

    xs = np.linspace(source_boundaries_pos_x[0], source_boundaries_pos_x[1], number_of_sources_x)
    ys = np.linspace(source_boundaries_pos_y[0], source_boundaries_pos_y[1], number_of_sources_y)

    xr2 = np.linspace(receivers_boundaries_pos_x[1][0], receivers_boundaries_pos_x[1][1], number_of_receivers_x)
    yr2 = np.linspace(receivers_boundaries_pos_y[1][0], receivers_boundaries_pos_y[1][1], number_of_receivers_y)

    receivers = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xr1 for y in yr1])
    receivers2 = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xr2 for y in yr2])

    receivers.append(receivers2)

   # receivers = receivers.getInsideDomain(box)
    shots = []

    for i in range(len(ys)):
        for j in range(len(xs)):

            srcpos = [xs[j], ys[i], source_depth]
            source = Source(srcpos, wavelet, dt)

            shot = Shot(source, receivers)
            shots.append(shot)

    return shots

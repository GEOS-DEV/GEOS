'''
Created on 9/02/2021

@author: macpro
'''
import numpy as np

import copy

from shot import *

from receiver import *
from source import *
import random
import math as m

#from matplotlib import pyplot as plt

'''acquisition sismique marine, comme schema PDF'''
def moving_acquisition(box, 
                       wavelet, 
                       nbsourcesx = 1, 
                       nbsourcesy = 1, 
                       nbreceiversx = 1, 
                       nbreceiversy = 1, 
                       lenRx = 0, 
                       lenRy = 0):
    
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
    ypos = np.linspace(ymin - lenRy/2,ymin + lenRy/2, nbreceiversy)
    
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
        
            # Define source location and type
            source = Source(box, srcpos, wavelet)
            source.setDomainFlag(box)
            
            receivers=receivers.getInsideDomain(box)
            
            if source.getDomainFlag()==1:
            # Create and store the shot
                shot = Shot(source, copy.deepcopy(receivers))
                shots.append(shot)

                
           #Uncomment this part to plot source and receivers  
            #plt.xlim(xmin,xmax)
            #plt.ylim(ymin,ymax)
            #plt.plot(source.getSourcePos()[0], source.getSourcePos()[1], 'bo')
            #for k in range(receivers.getNumberofReceivers()):
            #    plt.plot(receivers.getReceiverSetPos()[k][0], receivers.getReceiverSetPos()[k][1], 'ro')
            #plt.show()
    return shots

'''Create random acquisition'''
def random_acquisition(box, wavelet, nbsources, nbreceiversx, nbreceiversy):
    
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
        source = Source(box, srcpos, wavelet)
        source.setDomainFlag(box)
        
        if source.getDomainFlag()==1:
        # Create and store the shot
            shot = Shot(source, receivers)
            shots.append(shot)

    return shots


'''Create equiscaped acquisition (copy from Pysit)'''
def equispaced_acquisition(box, wavelet,
                           nbsourcesx,
                           nbsourcesy,
                           nbreceiversx,
                           nbreceiversy
                           ):


    xmin = box[0][0]
    xmax = box[0][1]
    
    ymin = box[1][0]
    ymax = box[1][1]
    
    zmin = box[2][0]
    zmax = box[2][1]


    xpos = np.linspace(xmin, xmax, nbreceiversx)
    ypos = np.linspace(ymin, ymax, nbreceiversy)
    
    receivers = ReceiverSet([Receiver([x, y, zmax]) for x in xpos for y in ypos])
    
    receivers=receivers.getInsideDomain(box)
    
    local_sources = nbsourcesx , nbsourcesy
    
    shots=list()

    for i in range(int(local_sources[0])):
        for j in range(int(local_sources[1])):

            idx = i + local_sources[0]
            jdx = j + local_sources[1]

            srcpos = [xmin + (xmax-xmin)*(idx+1.0)/(nbsourcesx+1.0), ymin + (ymax-ymin)*(jdx+1.0)/(nbsourcesy+1.0), zmax]
            
            # Define source location and type
            source = Source(box, srcpos, wavelet)
            source.setDomainFlag(box)
            
            if source.getDomainFlag()==1:
            # Create and store the shot
                shot = Shot(source, receivers)
                shots.append(shot)

    return shots

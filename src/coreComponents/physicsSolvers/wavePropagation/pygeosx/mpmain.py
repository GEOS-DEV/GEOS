'''
Created on 9/02/2021

@author: macpro
'''

from functions import *
from mesh import *
from acquisition import *
import pygeosx
import sys
from mpi4py import MPI
from segy import *
from print import *

import multiprocessing as mp


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():

    
    problem = pygeosx.initialize(0, sys.argv)
    # - This is a dummy initialization to be able to get the seismic acquisition/shot_list and split it between all proc
    # - Note that the current way to get a seismic acquisition will be replaced by the lecture of 
    #   a segy file ==> no longer need to do this dummy pygeosx.initialize()
    
    
    # Pseudo copy of GEOSX mesh
    mesh = initialize_pyMesh(problem) 
    
    maxT      = problem.get_wrapper("Events/maxTime").value()   
    dt        = calculDt(mesh)
    frequency = 10.0
    wavelet   = ricker(maxT[0], dt, frequency)
    
    #Set acquisition parameters 
    box             = mesh.getMinMaxBoundary()
    nb_source_x     = 5
    nb_source_y     = 2
    nb_receiver_x   = 3
    nb_receiver_y   = 4
    receiver_zone_x = 50
    receiver_zone_y = 90
    #Get seismic acquisition
    shot_list = moving_acquisition(box, wavelet, nb_source_x, nb_source_y, nb_receiver_x, nb_receiver_y, receiver_zone_x, receiver_zone_y) 
    
    
    """At this point we have defined our seismic acquisition/shot_list"""
    
    
    nb_proc = 2
    p = []
    nb_shot_m1 = len(shot_list)
    ind = 0
 
    #Loop over the process launch
    for i in range(nb_proc):
        nb_shot = int(nb_shot_m1/(nb_proc-i))
        
        #Add process with its dedicated shot_list
        p.append( mp.Process(target = shot_simul, 
                             args   = (shot_list[ind:ind + nb_shot], dt)) )
                             
        ind = ind + nb_shot      
        nb_shot_m1 = nb_shot_m1 - nb_shot
        
        #Start process
        p[i].start()
    
    
    for i in range(nb_proc):
        p[i].join()
    
    
if __name__ == "__main__":
    main()

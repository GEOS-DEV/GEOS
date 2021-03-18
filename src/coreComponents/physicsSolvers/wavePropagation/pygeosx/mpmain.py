'''
Created on 9/02/2021

@author: macpro
'''

from functions import *
import numpy as np
from mesh import *
from acquisition import *
import pygeosx
import sys
from mpi4py import MPI
from matplotlib import pyplot as plt
from segy import *
from print import *

import multiprocessing as mp


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def callback(matrix, array):
    pass



def main():

    #Dummy initialization to be able to get the shot list and split it between all proc
    problem = pygeosx.initialize(0, sys.argv)
    
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
    
    
    
    nb_proc = 2
   # mp.set_start_method('spawn')
    p=[]
    nb_shot_m1 = len(shot_list)
    ind = 0
 
    for i in range(nb_proc):
        nb_shot = int(nb_shot_m1/(nb_proc-i))
        
        p.append(mp.Process(target=shot_simul, args=(shot_list[ind:ind + nb_shot], dt) ))
        ind = ind + nb_shot
        
        nb_shot_m1 = nb_shot_m1 - nb_shot
        
        p[i].start()
    
    
    for i in range(nb_proc):
        p[i].join()
    
    
if __name__ == "__main__":
    main()

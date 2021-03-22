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


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():
    
    #Initialize GEOSX
    problem = pygeosx.initialize(rank, sys.argv)

    pygeosx.apply_initial_conditions()
    

#========================================================
#                 Initialization
#========================================================

    #Make a copy of mesh in Python object
    mesh = initialize_pyMesh(problem)      
    
    
    #Get Acoustic group
    acoustic_group  = problem.get_group("Solvers/acousticSolver")
    
    #Get Wrappers
    src_pos_geosx   = acoustic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    rcv_pos_geosx   = acoustic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)
    
    pressure_geosx  = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()
    
    pressure_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    pressure_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    pressure_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    
    outputSismoTrace = acoustic_group.get_wrapper("outputSismoTrace").value()
    
    dt_geosx        = problem.get_wrapper("Events/solverApplications/forceDt").value()
    maxT            = problem.get_wrapper("Events/maxTime").value()
    cycle_freq      = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle           = problem.get_wrapper("Events/python/lastCycle").value()
    cflFactor       = problem.get_wrapper("Solvers/acousticSolver/cflFactor").value()[0]
    curr_time       = problem.get_wrapper("Events/time").value()
    
    
    #Calculate dt using the inner sphere radius of the smallest element
    dt = calculDt(mesh)
    
    cycle_freq[0] = 1
    dt_geosx[0]   = dt       
    maxCycle      = int(maxT[0]/dt)
    maxT[0]       = (maxCycle+1) * dt 
    dt_cycle      = 1 
    
    

#========================================================
#                 Seismic Acquisition
#========================================================
    
    #Calculate the source function using the dt calculated previously
    frequency = 10.0
    wavelet   = ricker(maxT[0], dt, frequency)
    
    #Set acquisition parameters 
    box             = mesh.getMinMaxBoundary()
    nb_source_x     = 2
    nb_source_y     = 1
    nb_receiver_x   = 3
    nb_receiver_y   = 4
    receiver_zone_x = 50
    receiver_zone_y = 90
    #Get seismic acquisition
    shot_list = moving_acquisition(box, 
                                   wavelet, 
                                   nb_source_x, 
                                   nb_source_y, 
                                   nb_receiver_x, 
                                   nb_receiver_y, 
                                   receiver_zone_x, 
                                   receiver_zone_y) 
    
    nb_shot = len(shot_list)
    

#========================================================
#                      Solver 
#========================================================

    
    pressure_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers())) 
    
    ishot = 0
    print_shot_config(shot_list, ishot)
    
    #Set first source and receivers positions in GEOSX        
    src_pos          = shot_list[ishot].getSource().getCoord()       
    rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
    src_pos_geosx[0] = src_pos 
    rcv_pos_geosx.resize(len(rcv_pos_list))    
    rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]
    
    #Update shot flag
    shot_list[ishot].flagUpdate("In Progress")
    print_flag(shot_list)
    
    while (np.array([shot.getFlag() for shot in shot_list]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:
        #Save pressure
        if cycle[0] < (ishot+1) * maxCycle:
            pressure_at_receivers[cycle[0] - ishot * maxCycle, :] = pressure_geosx.to_numpy()[:] 
            
        else:
            pressure_at_receivers[maxCycle, :] = pressure_geosx.to_numpy()[:] 
            print_pressure(pressure_at_receivers, ishot)
           
            #Segy export and flag update
            if outputSismoTrace == 1 :
                export_to_segy(pressure_at_receivers, 
                               shot_list[0].getReceiverSet().getSetCoord(), 
                               ishot, 
                               dt_cycle)
           
            shot_list[ishot].flagUpdate("Done")
            
            #Reset time/pressure to 0
            curr_time[0] = 0.0
            pressure_nm1.to_numpy()[:] = 0.0
            pressure_n.to_numpy()[:]   = 0.0
            pressure_np1.to_numpy()[:] = 0.0
           
            #Increment shot
            ishot += 1             
            if ishot < nb_shot:
                print_shot_config(shot_list, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = shot_list[ishot].getSource().getCoord()       
                rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
                src_pos_geosx[0] = src_pos 
                rcv_pos_geosx.resize(len(rcv_pos_list))    
                rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]
                
                #Update shot flag
                shot_list[ishot].flagUpdate("In Progress") 
                
            print_flag(shot_list) 
                
                        
    
if __name__ == "__main__":
    main()

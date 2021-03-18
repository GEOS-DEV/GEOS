'''
Created on 8/02/2021

@author: macpro
'''

import math as m
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

    
    
def callback(matrix, array):
    pass
    
    
    
    
def shot_simul(shot_list, dt):
    problem = initialize()   
    
    do_shots(problem, shot_list, dt)
    
    
       

def initialize():
    problem = pygeosx.reinit(sys.argv)
    
    try:
        problem.get_group("Solvers/acousticSolver", None).register(callback)
    except AttributeError:
        pass

    pygeosx.apply_initial_conditions()
    
    return problem



def do_shots(problem, shot_list, dt):

    #Get Acoustic group
    acoustic_group  = problem.get_group("Solvers/acousticSolver")
    
    #Get Wrappers
    src_pos_geosx   = acoustic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    src_pos_geosx   = src_pos_geosx.to_numpy()
    
    rcv_pos_geosx_temp   = acoustic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx_temp.set_access_level(pygeosx.pylvarray.RESIZEABLE)
    
    pressure_geosx  = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()
    pressure_geosx  = pressure_geosx.to_numpy()
    
    pressure_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    pressure_nm1 = pressure_nm1.to_numpy()
    
    pressure_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    pressure_n = pressure_n.to_numpy()
    
    pressure_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)
    pressure_np1 = pressure_np1.to_numpy()
    
    dt_geosx        = problem.get_wrapper("Events/solverApplications/forceDt").value()
    maxT            = problem.get_wrapper("Events/maxTime").value()
    cycle_freq      = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle           = problem.get_wrapper("Events/python/lastCycle").value()
    cflFactor       = problem.get_wrapper("Solvers/acousticSolver/cflFactor").value()[0]
    curr_time       = problem.get_wrapper("Events/time").value()
    
    
    #Calculate dt using the inner sphere radius of the smallest element
    
    cycle_freq[0] = 1
    dt_geosx[0]   = dt       
    maxCycle      = int(maxT[0]/dt)
    maxT[0]       = (maxCycle+1) * dt 
    dt_cycle      = 1 
    
    pressure_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers())) 
    nb_shot = len(shot_list)
    
    ishot = 0
    print_shot_config(shot_list, ishot)
    
    #Set first source and receivers positions in GEOSX        
    src_pos          = shot_list[ishot].getSource().getCoord()       
    rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
    src_pos_geosx[0] = src_pos 
    rcv_pos_geosx_temp.resize(len(rcv_pos_list))    
    rcv_pos_geosx    = rcv_pos_geosx_temp.to_numpy()
    rcv_pos_geosx[:] = rcv_pos_list[:]
    
    #Update shot flag
    shot_list[ishot].flagUpdate("In Progress")
    print_flag(shot_list)
    
    while (np.array([shot.getFlag() for shot in shot_list]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:
        #Extract pressure
        pressure_geosx = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()
        pressure_geosx = pressure_geosx.to_numpy()
        
        #Save pressure
        if cycle[0] < (ishot+1) * maxCycle:
            pressure_at_receivers[cycle[0] - ishot * maxCycle, :] = pressure_geosx[:] 
            
        else:
            pressure_at_receivers[maxCycle, :] = pressure_geosx[:] 
            print_pressure(pressure_at_receivers, ishot)
           
            #Segy export and flag update
            export_to_segy(pressure_at_receivers, shot_list[0].getReceiverSet().getSetCoord(), ishot, dt_cycle)
            shot_list[ishot].flagUpdate("Done")
            
            #Reset time/pressure to 0
            curr_time[0]    = 0.0
            pressure_nm1[:] = 0.0
            pressure_n[:]   = 0.0
            pressure_np1[:] = 0.0
           
            #Increment shot
            ishot += 1             
            if ishot < nb_shot:
                print_shot_config(shot_list, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = shot_list[ishot].getSource().getCoord()       
                rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord() 
    
                src_pos_geosx[0] = src_pos 
                rcv_pos_geosx_temp.resize(len(rcv_pos_list))    
                rcv_pos_geosx    = rcv_pos_geosx_temp.to_numpy()
                rcv_pos_geosx[:] = rcv_pos_list[:]
                
                #Update shot flag
                shot_list[ishot].flagUpdate("In Progress") 
                
            print_flag(shot_list) 




def ricker(maxT, dt, f0):
    
    T0 = 1.0/f0;
    fi = np.zeros(int(maxT/dt))
    
    for t in range(int(maxT/dt)):
        t0 = dt*t
        if t0 <= -0.9*T0 or t0 >= 2.9*T0:
            fi[t] = 0.0;
        else:
            tmp      = f0*t0-1.0
            f0tm1_2  = 2*(tmp*np.pi)*(tmp*np.pi)
            gaussian = m.exp( -f0tm1_2 )
            fi[t]    = -(t0-1)*gaussian
               
    return fi



#Calculate dt using order of space discretization method, Wave velocity, and the radius of the included sphere in element 
def calculDt(mesh):
   
    if mesh.getOrd()==1:
        nx = 1
        ny = 2
        nz = 4
    elif mesh.getOrd()==3:
        nx = 3
        ny = 12
        nz = 48
    elif mesh.getOrd()==5:
        nx = 5
        ny = 30
        nz = 180
    

    h    = np.linalg.norm(mesh.getElem_List()[0].getNode_List()[0].getCoords() - mesh.getElem_List()[0].getNode_List()[nx].getCoords())/2
    Vmax = mesh.getElem_List()[0].getSpeed()
    V    = mesh.getElem_List()[0].getVolume()
    
    for elem in (mesh.getElem_List()):
        
        if elem.getVolume()<V: 
        
            xradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nx].getCoords())/2
            yradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[ny].getCoords())/2
            zradius = np.linalg.norm(elem.getNode_List()[0].getCoords() - elem.getNode_List()[nz].getCoords())/2
            
            if xradius < h:
                h = xradius
            if yradius < h:
                h = yradius
            if zradius < h:
                h = zradius
            
        if Vmax < elem.getSpeed():
            Vmax = elem.getSpeed()
            
        
    dt = h/(Vmax*mesh.ord)           
            
    return dt
    

            

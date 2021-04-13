'''
Created on 8/02/2021

@author: macpro
'''

import numpy as np
import pygeosx

from segyManager import *
from print import *


def shot_simul(rank, xml, shot_list, tracePath):
    """
    Parameters
    ----------

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    dt : float
        Time step for the solver
    """

    problem = initialize(rank, xml)
    do_shots(rank, problem, shot_list, tracePath)



def initialize(rank, xml):
    """ Grouping of pygeox initializations

    Return
    ------
    problem :
        The pygeosx ProblemManager

    Notes
    -----
    Need to give MPI rank at this point for initialization.
    Conflict with first initialization to get the list of shots
    """

    problem = pygeosx.initialize(rank, xml)
    pygeosx.apply_initial_conditions()

    return problem



def do_shots(rank, problem, shot_list, tracePath):
    """ Given a GEOSX problem, a list of shots, and a time step,
        solve wave eqn with different configurations

    Parameters
    ----------
    problem :
        The pygeosx ProblemManager

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    dt : float
        Time step for the solver

    Notes
    -----
    Export pressure to segy depending on output flag defined in XML
    """

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

    dt = shot_list[0].getSource().getTimeStep()
    cycle_freq[0] = 1
    dt_geosx[0]   = dt
    maxCycle      = int(maxT[0]/dt)
    maxT[0]       = (maxCycle+1) * dt
    dt_cycle      = 1 #to be changed

    pressure_at_receivers = np.zeros((maxCycle+1, shot_list[0].getReceiverSet().getNumberOfReceivers()))
    nb_shot = len(shot_list)

    ishot = 0
    if rank:
        print_shot_config(shot_list, ishot)

    #Set first source and receivers positions in GEOSX
    src_pos          = shot_list[ishot].getSource().getCoord()
    rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord()

    src_pos_geosx.to_numpy()[0] = src_pos
    rcv_pos_geosx.resize(len(rcv_pos_list))
    rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

    #Update shot flag
    shot_list[ishot].flagUpdate("In Progress")
    if rank:
        print_flag(shot_list)

    while (np.array([shot.getFlag() for shot in shot_list]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:
        #Save pressure
        if cycle[0] < (ishot+1) * maxCycle:
            pressure_at_receivers[cycle[0] - ishot * maxCycle, :] = pressure_geosx.to_numpy()[:]

        else:
            pressure_at_receivers[maxCycle, :] = pressure_geosx.to_numpy()[:]
            if rank:
                print_pressure(pressure_at_receivers, ishot)

            #Segy export and flag update
            if rank:
                if outputSismoTrace == 0 :
                    export_to_segy(pressure_at_receivers,
                                   shot_list[ishot].getSource().getCoord(),
                                   shot_list[ishot].getReceiverSet().getSetCoord(),
                                   ishot,
                                   dt_cycle,
                                   tracePath)

            shot_list[ishot].flagUpdate("Done")

            #Reset time/pressure to 0
            curr_time[0]    = 0.0
            pressure_nm1.to_numpy()[:] = 0.0
            pressure_n.to_numpy()[:]   = 0.0
            pressure_np1.to_numpy()[:] = 0.0

            #Increment shot
            ishot += 1
            if ishot < nb_shot:
                if rank:
                    print_shot_config(shot_list, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = shot_list[ishot].getSource().getCoord()
                rcv_pos_list     = shot_list[ishot].getReceiverSet().getSetCoord()

                src_pos_geosx.to_numpy()[0] = src_pos
                rcv_pos_geosx.resize(len(rcv_pos_list))
                rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

                #Update shot flag
                shot_list[ishot].flagUpdate("In Progress")

            if rank:
                print_flag(shot_list)

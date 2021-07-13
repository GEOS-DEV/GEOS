'''
Created on 8/02/2021

@author: macpro
'''

import numpy as np
import pygeosx
from math import ceil

from segyManager import *
from fileManager import rootPath
import os
from print import *


def acoustic_shots(rank, problem, acquisition, tracePath):
    """ Given a GEOSX problem, a list of shots, and a time step,
        solve wave eqn with different configurations

    Parameters
    ----------
    rank : int
        MPI rank

    problem :
        The pygeosx ProblemManager

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    tracePath : string
        Directory where to output .sgy files

    Notes
    -----
    Export pressure to segy depending on output flag defined in XML
    """

    #Get Acoustic group
    acoustic_group  = problem.get_group("Solvers/acousticSolver")


    #Get Wrappers
    src_pos_geosx  = acoustic_group.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    rcv_pos_geosx  = acoustic_group.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)

    pressure_geosx = acoustic_group.get_wrapper("pressureNp1AtReceivers").value()

    pressure_nm1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_n = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_np1 = problem.get_wrapper("domain/MeshBodies/mesh/Level0/nodeManager/pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    nsamples = acoustic_group.get_wrapper("nSampleSismoTrace").value()
    outputSismoTrace = acoustic_group.get_wrapper("outputSismoTrace").value()
    indexSismoTrace = acoustic_group.get_wrapper("indexSismoTrace").value()

    dt_geosx   = problem.get_wrapper("Events/solverApplications/forceDt").value()
    dt_sismo   = acoustic_group.get_wrapper("dtSismoTrace").value()
    maxT       = problem.get_wrapper("Events/maxTime").value()
    cycle_freq = problem.get_wrapper("Events/python/cycleFrequency").value()
    cycle      = problem.get_wrapper("Events/python/lastCycle").value()
    cycle2     = problem.get_wrapper("Events/cycle").value()
    curr_time  = problem.get_wrapper("Events/time").value()

    dt            = acquisition.dt
    dt_geosx[0]   = dt
    maxCycle      = int(ceil(maxT[0]/dt))
    cycle_freq[0] = maxCycle
    maxT[0]       = (maxCycle+1)*dt
    nsamples[0]   = int(maxT[0]/dt_sismo)
    
    pressure_at_receivers = np.zeros((nsamples[0], acquisition.shot[0].receivers.n))
    nb_shot = len(acquisition.shot)

    if outputSismoTrace == 1 :
        if rank==0:
            create_segy(acquisition.shot, "pressure", nsamples[0], acquisition.output)


    ishot = 0
    if rank==0:
        print_shot_config(acquisition.shot, ishot)

    #Set first source and receivers positions in GEOSX
    src_pos      = acquisition.shot[ishot].source.coords
    rcv_pos_list = [receiver.coords for receiver in acquisition.shot[ishot].receivers.receiver_list]

    src_pos_geosx.to_numpy()[0] = src_pos
    rcv_pos_geosx.resize(len(rcv_pos_list))
    rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

    #Update shot flag
    acquisition.shot[ishot].flag = "In Progress"
    if rank==0:
        print_flag(acquisition.shot)

    while (np.array([shot.flag for shot in acquisition.shot]) == "Done").all() != True and pygeosx.run() != pygeosx.COMPLETED:
        #Save pressure
        if cycle[0] !=0 :
            pressure_at_receivers[:, :] = pressure_geosx.to_numpy().transpose()

            #Segy export and flag update
            if outputSismoTrace == 1 :
                segyFile = os.path.join(acquisition.output, "pressure_Shot"+ acquisition.shot[ishot].id + ".sgy")
                export_to_segy(pressure_at_receivers,
                               acquisition.shot[ishot].receivers.receiver_list,
                               segyFile)

            acquisition.shot[ishot].flag = "Done"

            #Reset time to -dt and pressure to 0
            curr_time[0]               = -dt
            cycle_freq[0]              = maxCycle + 1
            #cycle2[0]                  = 0
            pressure_nm1.to_numpy()[:] = 0.0
            pressure_n.to_numpy()[:]   = 0.0
            pressure_np1.to_numpy()[:] = 0.0
            indexSismoTrace[0]         = 0

            #Increment shot
            ishot += 1
            if ishot < nb_shot:
                if rank==0:
                    print_shot_config(acquisition.shot, ishot)

                #Set new receivers and source positions in GEOSX
                src_pos          = acquisition.shot[ishot].source.coords
                rcv_pos_list     = [receiver.coords for receiver in acquisition.shot[ishot].receivers.receiver_list]

                src_pos_geosx.to_numpy()[0] = src_pos
                rcv_pos_geosx.resize(len(rcv_pos_list))
                rcv_pos_geosx.to_numpy()[:] = rcv_pos_list[:]

                #Update shot flag
                acquisition.shot[ishot].flag = "In Progress"

            if rank==0:
                print_flag(acquisition.shot)

import pygeosx

import sys
import os
import numpy as np
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.segy import exportToSegy
from seismicUtilities.acoustic import updateSourceAndReceivers, setTimeVariables, print_group
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():
    #Set acquisition

    acquisition = EQUISPACEDAcquisition(boundary=[[0,2000],[0,2000],[0,2000]],
                                        dt=0.001,
                                        velocity_model=1500,
                                        start_source_pos    = [1001, 1001],
                                        end_source_pos      = [1001, 1001],
                                        start_receivers_pos = [[21, 1001]],
                                        end_receivers_pos   = [[1951, 1001]],
                                        number_of_sources   = 1,
                                        number_of_receivers = 150,
                                        source_depth = 99,
                                        receivers_depth = 49)

    acquisition.add_xml(sys.argv[2])
    #acquisition.limitedAperture(500)
    #acquisition.calculDt()


    maxTime = 2.0
    nbSeismo = 500
    outputWaveFieldInterval = 10


    result = acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm)



def acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm):
    #Loop over the shots
    segyList = []
    ishot=0
    rank = comm.Get_rank()
#========================================================
    #INITIALIZATION
#========================================================
    #Loop over the shots
    for shot in acquisition.shots:
        sys.argv[2] = shot.xml
        dt = shot.dt
        dtSeismoTrace = maxTime/(nbSeismo - 1)

        #Initialize GEOSX
        if ishot == 0:
            problem = pygeosx.initialize(rank, sys.argv)
        else:
            problem = pygeosx.reinit(sys.argv)
        #Get solver, get hdf5 collection/outputs wrappers
        acousticSolver = problem.get_group("/Solvers/acousticSolver")

        #Set times variables
        setTimeVariables(problem, maxTime, dt, dtSeismoTrace)

#===================================================
        #FORWARD
#===================================================

        #Get view on pressure at receivers locations
        pressureAtReceivers = acousticSolver.get_wrapper("pressureNp1AtReceivers").value()

        #update source and receivers positions
        updateSourceAndReceivers(acousticSolver, shot.sources.source_list, shot.receivers.receivers_list)
        pygeosx.apply_initial_conditions()

        shot.flag = "In Progress"

        time = 0.0
        i = 0

        if rank == 0 :
            print("\nForward")

        #Time loop
        while time < maxTime:
            if rank == 0:
                print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", i+1)

            #Execute one time step
            acousticSolver.execute(time, dt)

            time += dt
            i += 1

        #Export pressure at receivers to segy
        pressure = np.array(pressureAtReceivers.to_numpy())
        exportToSegy(table = pressure,
                     shot = shot,
                     filename = "seismo",
                     directory = "dataTest",
                     rank = rank)

        pygeosx._finalize()

if __name__ == "__main__":
    main()

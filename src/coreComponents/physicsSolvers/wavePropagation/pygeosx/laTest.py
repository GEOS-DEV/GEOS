import pygeosx
import sys
from acquisition import EQUISPACEDAcquisition
from print import *
import segyio
import numpy as np
from segyManager import *
from mpi4py import MPI

import h5py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def recomputeSourceAndReceivers(solver, source, receivers):
    updateSourceAndReceivers(solver, source, receivers)

    solver.initPostInitialConditions()


def updateSourceAndReceivers(solver, source, receivers):
    src_pos_geosx = solver.get_wrapper("sourceCoordinates").value()
    src_pos_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    rcv_pos_geosx = solver.get_wrapper("receiverCoordinates").value()
    rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)


    src_pos_geosx.to_numpy()[0] = source.coords
    rcv_pos_geosx.resize(receivers.n)
    rcv_pos = [receiver.coords for receiver in receivers.receivers_list]
    rcv_pos_geosx.to_numpy()[:] = rcv_pos[:]

    solver.postProcessInput()



def resetWaveField(group):
    group.get_wrapper("Solvers/acousticSolver/indexSeismoTrace").value()[0] = 0
    nodeManagerPath = "domain/MeshBodies/mesh/Level0/nodeManager/"

    pressure_nm1 = group.get_wrapper(nodeManagerPath + "pressure_nm1").value()
    pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_n = group.get_wrapper(nodeManagerPath + "pressure_n").value()
    pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_np1 = group.get_wrapper(nodeManagerPath + "pressure_np1").value()
    pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_geosx = group.get_wrapper("Solvers/acousticSolver/pressureNp1AtReceivers").value()
    pressure_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

    pressure_nm1.to_numpy()[:] = 0.0
    pressure_n.to_numpy()[:]   = 0.0
    pressure_np1.to_numpy()[:] = 0.0
    pressure_geosx.to_numpy()[:] = 0.0


def setTimeVariables(problem, maxTime, dtSeismoTrace):
    problem.get_wrapper("Events/maxTime").value()[0] = maxTime
    problem.get_wrapper("/Solvers/acousticSolver/dtSeismoTrace").value()[0] = dtSeismoTrace





def main():
    acquisition = EQUISPACEDAcquisition(boundary=[[0,2000],[0,2000],[0,2000]],
                                        dt=0.005,
                                        velocity_model=1500,
                                        start_source_pos    = [501, 1001],
                                        end_source_pos      = [1501, 1001],
                                        start_receivers_pos = [[21, 1001]],
                                        end_receivers_pos   = [[1951, 1001]],
                                        number_of_sources   = 1,
                                        number_of_receivers = 50,
                                        source_depth = 1901,
                                        receivers_depth = 1951)

    acquisition.add_xml(sys.argv[2])
    acquisition.limitedAperture(500)
    comm.Barrier()
    acquisition.calculDt()



    maxTime = 2.0
    outputInterval = 5

    #Loop over the shots
    ishot=0
    for shot in acquisition.shots:

        dt = shot.dt
        dtSeismoTrace = dt * outputInterval
        nsamples = int(maxTime/(dt*outputInterval))+1

        if rank == 0:
            create_segy(shot, "pressure", nsamples, acquisition.output)

        #Initialize GEOSX, get Acoustic Solver, get HDF5 Outputs
        if ishot == 0:
            problem = pygeosx.initialize(rank, sys.argv)
        else:
            problem = pygeosx.reinit(sys.argv)

        acousticSolver = pygeosx.pysolver.Solver("/Solvers/acousticSolver")
        hdf5 = pygeosx.pyhdf5.HDF5()

        #Get view on pressure at receivers locations
        pressure_geosx = acousticSolver.get_wrapper("pressureNp1AtReceivers").value()

        #Set time variables and update source and receivers positions
        setTimeVariables(problem, maxTime, dtSeismoTrace)
        updateSourceAndReceivers(acousticSolver, shot.source, shot.receivers)
        pygeosx.apply_initial_conditions()

        shot.flag = "In Progress"

        time = 0.0
        i = 0
        while time < maxTime:
            if rank == 0:
                print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", i+1)

            acousticSolver.explicitStep(time, dt)

            time += dt
            i += 1
            #Collect waveField values
            if i % 100 == 0:
                hdf5.collect("waveField", time, dt)

        #Export pressure at receivers positions to .segy file
        segyFile = os.path.join(acquisition.output, "pressure_Shot"+ shot.id + ".sgy")
        export_to_segy(pressure_geosx.to_numpy(),
                       shot.receivers.receivers_list,
                       segyFile)

        #Output waveField values
        hdf5.output("waveField", time, dt)

        shot.flag = "Done"
        resetWaveField(problem)

        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1

    #Uncomment to check if hdf5 file is correctly filled
    '''
    with h5py.File("waveField.hdf5", "r") as f:
    # List all groups
        print("Keys: %s" % f.keys())
        a_group_key = list(f.keys())[0]

        # Get the data
        data = list(f[a_group_key])

    print(len(data[0])*len(data))
    '''

    pygeosx._finalize()



if __name__ == "__main__":
    main()

import pygeosx

import sys
import numpy as np
import h5py
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.segy import exportToSegy
from seismicUtilities.acoustic import updateSourceAndReceivers, updateSourceValue, residualLinearInterpolation, resetWaveField, setTimeVariables
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():
    #Set acquisition

    acquisition = EQUISPACEDAcquisition(boundary=[[0,2000],[0,2000],[0,2000]],
                                        dt=0.005,
                                        velocity_model=1500,
                                        start_source_pos    = [501, 1001],
                                        end_source_pos      = [1501, 1001],
                                        start_receivers_pos = [[21, 1001]],
                                        end_receivers_pos   = [[1951, 1001]],
                                        number_of_sources   = 1,
                                        number_of_receivers = 2,
                                        source_depth = 1901,
                                        receivers_depth = 1951)
    """
    acquisition = EQUISPACEDAcquisition(boundary=[[0,13520],[0,13520],[0,4200]],
                                dt=0.002,
                                velocity_model="/home/julien/codes/SEP_REDUCE_Model/338x338x105_velModel.geos",
                                start_source_pos    = [7001, 7001],
                                end_source_pos      = [12001, 7001],
                                start_receivers_pos = [[21, 7001]],
                                end_receivers_pos   = [[13501, 7001]],
                                number_of_sources   = 1,
                                number_of_receivers = 676,
                                source_depth = 4099,
                                receivers_depth = 4149)
    """
    acquisition.add_xml(sys.argv[2])
    #acquisition.limitedAperture(500)
    #acquisition.calculDt()


    maxTime = 0.5
    nbSeismo = 10
    dtSeismoTrace = maxTime/nbSeismo
    outputWaveFieldInterval = 30


    result = acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, rank)



def acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, rank=0):
    #Loop over the shots
    segyList = []
    ishot=0

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
        collection = problem.get_group("/Tasks/waveFieldCollection")
        output = problem.get_group("Outputs/waveFieldOutput")

        #Set times variables
        setTimeVariables(problem, maxTime, dt, dtSeismoTrace)

#===================================================
        #FORWARD
#===================================================
        #Update hdf5 output filename
        output.setOutputName("forwardWaveField"+shot.id)

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

            #Collect/export waveField values to hdf5
            if i % outputWaveFieldInterval == 0:
                collection.collect(time, dt)
                output.output(time, dt)

        #Export pressure at receivers to segy
        pressure = np.array(pressureAtReceivers.to_numpy())
        exportToSegy(table = pressure,
                     shot = shot,
                     filename = "pressure_Shot"+shot.id,
                     directory = acquisition.output,
                     rank = rank)


        #Residual computation
        residual = residualLinearInterpolation(pressure, maxTime, dt, dtSeismoTrace)

        #Reset waveField
        resetWaveField(problem)

#==================================================
        #BACKWARD
#==================================================
        #initialize adjoint problem
        updateSourceAndReceivers(acousticSolver, shot.receivers.receivers_list)
        updateSourceValue(acousticSolver, residual)

        #Update hdf5 output filename
        output.setOutputName("backwardWaveField"+shot.id)
        output.reinit()

        if rank == 0:
            print("\nbackward")

        #Reverse time loop
        while time > 0:
            if rank == 0:
                print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", i)

            #Execute one time step backward
            acousticSolver.execute(time, dt)

            #Collect/export waveField values to hdf5
            if i % outputWaveFieldInterval == 0:
                collection.collect(time, dt)
                output.output(time, dt)

            time -= dt
            i -= 1


#==================================================
        #PARTIAL GRADIENT
#==================================================

        if rank == 0:
            h5f = h5py.File("forwardWaveField"+shot.id+".hdf5", "r")
            h5b = h5py.File("backwardWaveField"+shot.id+".hdf5", "r")
            keys = list(h5f.keys())

            h5p = h5py.File("partialGradient"+shot.id+".hdf5", "w")
            h5p.create_dataset("partialGradient_np1", data = h5f[keys[0]], dtype='d')
            h5p.create_dataset("partialGradient_np1 ReferencePosition", data = h5f[keys[1]])
            h5p.create_dataset("partialGradient_np1 Time", data = h5f[keys[2]])
            keyp = list(h5p.keys())[0]

            n = len( list( h5f[keys[0]] ) )
            for i in range(n):
                # Get the data
                h5p[keyp][i][:] *= h5b[keys[0]][n-1-i][:]


        #Update shot flag
        shot.flag = "Done"
        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1

    pygeosx._finalize()

if __name__ == "__main__":
    main()

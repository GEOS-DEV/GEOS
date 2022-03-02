import pygeosx

import sys
import os
import numpy as np
import h5py
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.segy import exportToSegy
from seismicUtilities.acoustic import updateSourceAndReceivers, updateSourceValue, residualLinearInterpolation, resetWaveField, setTimeVariables, computeFullGradient, computeResidual
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


    maxTime = 2.0
    nbSeismo = 500
    outputWaveFieldInterval = 50


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
        collectionNp1 = problem.get_group("/Tasks/waveFieldNp1Collection")
        collectionN = problem.get_group("/Tasks/waveFieldNCollection")
        collectionNm1 = problem.get_group("/Tasks/waveFieldNm1Collection")
        outputNp1 = problem.get_group("Outputs/waveFieldNp1Output")
        outputN = problem.get_group("Outputs/waveFieldNOutput")
        outputNm1 = problem.get_group("Outputs/waveFieldNm1Output")
        #Set times variables
        setTimeVariables(problem, maxTime, dt, dtSeismoTrace)

#===================================================
        #FORWARD
#===================================================
        #Update hdf5 output filename
        directory = "partialGradient"
        if rank == 0:
            if os.path.exists(directory):
                pass
            else:
                os.mkdir(directory)

        outputNp1.setOutputName(os.path.join(directory,"forwardWaveFieldNp1"+shot.id))
        outputN.setOutputName(os.path.join(directory,"forwardWaveFieldN"+shot.id))
        outputNm1.setOutputName(os.path.join(directory,"forwardWaveFieldNm1"+shot.id))

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
                collectionNp1.collect(time, dt)
                outputNp1.output(time, dt)

                collectionN.collect(time, dt)
                outputN.output(time, dt)

                collectionNm1.collect(time, dt)
                outputNm1.output(time, dt)

        #Export pressure at receivers to segy
        pressure = np.array(pressureAtReceivers.to_numpy())
        exportToSegy(table = pressure,
                     shot = shot,
                     filename = "seismo",
                     directory = acquisition.output,
                     rank = rank)


        #Residual computation
        residualTemp = computeResidual("dataTest/seismo_Shot"+shot.id+".sgy", pressure)
        residual = residualLinearInterpolation(residualTemp, maxTime, dt, dtSeismoTrace)

        #Reset waveField
        resetWaveField(problem)

#==================================================
        #BACKWARD
#==================================================
        #initialize adjoint problem
        updateSourceAndReceivers(acousticSolver, shot.receivers.receivers_list)
        updateSourceValue(acousticSolver, residual)

        #Update hdf5 output filename
        outputNp1.setOutputName(os.path.join(directory,"backwardWaveFieldNp1"+shot.id))
        outputNp1.reinit()

        outputN.setOutputName(os.path.join(directory,"backwardWaveFieldN"+shot.id))
        outputN.reinit()

        outputNm1.setOutputName(os.path.join(directory,"backwardWaveFieldNm1"+shot.id))
        outputNm1.reinit()

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
                collectionNp1.collect(time, dt)
                outputNp1.output(time, dt)

                collectionN.collect(time, dt)
                outputN.output(time, dt)

                collectionNm1.collect(time, dt)
                outputNm1.output(time, dt)

            time -= dt
            i -= 1


#==================================================
        #PARTIAL GRADIENT
#==================================================
        if rank == 0:
            if os.path.exists(directory):
                pass
            else:
                os.mkdir(directory)

            h5fNp1 = h5py.File(os.path.join(directory,"forwardWaveFieldNp1"+shot.id+".hdf5"), "r+")
            h5bNp1 = h5py.File(os.path.join(directory,"backwardWaveFieldNp1"+shot.id+".hdf5"), "r+")

            h5fN = h5py.File(os.path.join(directory,"forwardWaveFieldN"+shot.id+".hdf5"), "r+")
            h5bN = h5py.File(os.path.join(directory,"backwardWaveFieldN"+shot.id+".hdf5"), "r+")

            h5fNm1 = h5py.File(os.path.join(directory,"forwardWaveFieldNm1"+shot.id+".hdf5"), "r+")
            h5bNm1 = h5py.File(os.path.join(directory,"backwardWaveFieldNm1"+shot.id+".hdf5"), "r+")

            keysNp1 = list(h5fNp1.keys())
            keysN = list(h5fN.keys())
            keysNm1 = list(h5fNm1.keys())

            for i in range(len(keysNp1)-2):
                if ("ReferencePosition" in keysNp1[i]) or ("Time" in keysNp1[i]):
                    del h5fNp1[keysNp1[i]]
                    del h5bNp1[keysNp1[i]]

                    del h5fN[keysNp1[i]]
                    del h5bN[keysNp1[i]]

                    del h5fNm1[keysNp1[i]]
                    del h5bNm1[keysNp1[i]]

            keysNp1 = list(h5fNp1.keys())

            h5p = h5py.File(os.path.join(directory,"partialGradient"+shot.id+".hdf5"), "w")

            h5p.create_dataset("ReferencePosition", data = h5fNp1[keysNp1[-2]][0], dtype='d')
            h5p.create_dataset("Time", data = h5fNp1[keysNp1[-1]], dtype='d')
            h5p.create_dataset("partialGradient", data = np.zeros(h5bNp1[keysNp1[0]][0].shape[0]), dtype='d')


            keyp = list(h5p.keys())

            n = len( list( h5fNp1[keysNp1[0]] ) )
            for i in range(n):
                #h5p["partialGradient"][-i-1, :] *= -(h5f["pressure_np1"][i, :] - 2*h5f["pressure_n"][i, :] + h5f["pressure_nm1"][i, :]) / (dt*dt)
                h5p["partialGradient"][:] -= (h5fNp1["pressure_np1"][i, :] - 2*h5fN["pressure_n"][i, :] + h5fNm1["pressure_nm1"][i, :]) * h5bN["pressure_n"][-i-1, :] / (dt*dt)


            h5fNp1.close()
            h5bNp1.close()
            h5fN.close()
            h5bN.close()
            h5fNm1.close()
            h5bNm1.close()
            h5p.close()

            os.remove(os.path.join(directory,"forwardWaveFieldNp1"+shot.id+".hdf5"))
            os.remove(os.path.join(directory,"backwardWaveFieldNp1"+shot.id+".hdf5"))
            os.remove(os.path.join(directory,"forwardWaveFieldN"+shot.id+".hdf5"))
            os.remove(os.path.join(directory,"backwardWaveFieldN"+shot.id+".hdf5"))
            os.remove(os.path.join(directory,"forwardWaveFieldNm1"+shot.id+".hdf5"))
            os.remove(os.path.join(directory,"backwardWaveFieldNm1"+shot.id+".hdf5"))

        #Update shot flag
        shot.flag = "Done"
        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1
        comm.Barrier()

    if rank==0:
        computeFullGradient(directory, acquisition)
    pygeosx._finalize()

if __name__ == "__main__":
    main()

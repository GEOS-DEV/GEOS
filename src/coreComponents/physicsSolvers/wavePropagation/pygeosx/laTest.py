import sys
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.fwi import forward,                       \
                                 backward,                      \
                                 computeFullGradient,           \
                                 computePartialGradient,        \
                                 computePartialCostFunction,    \
                                 computeFullCostFunction,       \
                                 computeResidual
from seismicUtilities.AcousticSolver import AcousticSolver

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
    #acquisition.limitedAperture(2000)
    #acquisition.calculDt()

    maxTime = 2.0
    nbSeismo = 500
    outputWaveFieldInterval = 50


    result = acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm)

def acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm):
    #Loop over the shots
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

        if ishot == 0:
            acousticSolver = AcousticSolver(sys.argv[2],
                                            dt,
                                            maxTime,
                                            dtSeismoTrace)

            acousticSolver.initialize(rank)

        else:
            acousticSolver = AcousticSolver(sys.argv[2],
                                            dt,
                                            maxTime,
                                            dtSeismoTrace)

            acousticSolver.reinitialize()

        gradDir = "partialGradient"
        acousticSolver.updateOutputsName(directory=gradDir,
                                         filenames=["forwardWaveFieldNp1_"+shot.id,
                                                    "forwardWaveFieldN_"+shot.id,
                                                    "forwardWaveFieldNm1_"+shot.id],)


        acousticSolver.apply_initial_conditions()

        acousticSolver.updateSourceAndReceivers(shot.sources.source_list, shot.receivers.receivers_list)
        #Get view on pressure at receivers locations
        pressureAtReceivers = acousticSolver.getPressureAtReceivers()

#===================================================
        #FORWARD
#===================================================
        residual = forward(acousticSolver, shot, outputWaveFieldInterval, residual_flag=True, datafile="dataTest/seismo_Shot"+shot.id+".sgy", rank=rank)

#==================================================
        #BACKWARD
#==================================================
        acousticSolver.updateSourceAndReceivers(sources_list = shot.receivers.receivers_list)
        acousticSolver.updateSourceValue(residual)
        acousticSolver.updateOutputsName(directory=gradDir,
                                         filenames=["backwardWaveFieldNp1_"+shot.id,
                                                    "backwardWaveFieldN_"+shot.id,
                                                    "backwardWaveFieldNm1_"+shot.id],
                                         reinit=True)

        backward(acousticSolver, shot, outputWaveFieldInterval, rank)

#==================================================
        #PARTIAL GRADIENT
#==================================================
        if rank == 0:
            computePartialGradient(gradDir, shot)

        #Update shot flag
        shot.flag = "Done"
        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1
        comm.Barrier()

#==================================================
        #FULL GRADIENT
#==================================================
    costDir = "partialCostFunction"
    if rank==0:
        fullCostFunction = computeFullCostFunction(costDir, acquisition)
        computeFullGradient(gradDir, acquisition)

    acousticSolver.finalize()

if __name__ == "__main__":
    main()

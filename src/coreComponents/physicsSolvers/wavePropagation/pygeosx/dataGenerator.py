import sys
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.fwi import forward,                       \
                                 backward,                      \
                                 computeFullGradient,           \
                                 computePartialGradient,        \
                                 computePartialCostFunction,    \
                                 computeFullCostFunction,       \
                                 computeResidual
from seismicUtilities.segy import exportToSegy
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
    #acquisition.limitedAperture(500)
    #acquisition.calculDt()

    maxTime = 2.0
    nbSeismo = 500
    outputWaveFieldInterval = 50


    acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm)

def acousticShot(maxTime, nbSeismo, outputWaveFieldInterval, acquisition, comm):
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

        acousticSolver = AcousticSolver(sys.argv[2],
                                        dt,
                                        maxTime,
                                        dtSeismoTrace)

        acousticSolver.initialize(rank)
        acousticSolver.apply_initial_conditions()

        acousticSolver.updateSourceAndReceivers(shot.sources.source_list, shot.receivers.receivers_list)
#===================================================
        #FORWARD
#===================================================

        residual = forward(acousticSolver, shot, outputWaveFieldInterval, rank)

        pressure = acousticSolver.getPressureAtReceivers()
        exportToSegy(table = pressure,
                     shot = shot,
                     filename = "seismo",
                     directory = "dataTest",
                     rank = rank)

    acousticSolver.finalize()

if __name__ == "__main__":
    main()

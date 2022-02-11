import pygeosx
import sys
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.segy import exportToSegy
from seismicUtilities.acoustic import updateSourceAndReceivers, recomputeSourceAndReceivers, resetWaveField, setTimeVariables
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
                                        number_of_sources   = 2,
                                        number_of_receivers = 50,
                                        source_depth = 1901,
                                        receivers_depth = 1951)

    acquisition.add_xml(sys.argv[2])
    #acquisition.calculDt()

    #Initialize GEOSX, get Acoustic Solver, get HDF5 Outputs
    problem = pygeosx.initialize(rank, sys.argv)
    acousticSolver = problem.get_group("/Solvers/acousticSolver")

    #Set time variables
    maxTime = 2.0
    outputInterval = 10
    dt = acquisition.shots[0].dt
    dtSeismoTrace = dt * outputInterval
    nsamples = int(maxTime/(dt*outputInterval)) +1

    setTimeVariables(problem, maxTime, dtSeismoTrace)

    #Get view on pressure at receivers locations
    pressureAtReceivers = acousticSolver.get_wrapper("pressureNp1AtReceivers").value()

    #Update source and receivers positions
    updateSourceAndReceivers(acousticSolver, acquisition.shots[0].source, acquisition.shots[0].receivers)
    pygeosx.apply_initial_conditions()

    #Loop over the shots
    ishot=0
    for shot in acquisition.shots:
        if ishot > 0:
            recomputeSourceAndReceivers(acousticSolver, shot.source, shot.receivers)

        shot.flag = "In Progress"

        time = 0.0
        i = 0

        #Time loop
        while time < maxTime:
            if rank == 0:
                print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", i+1)

            acousticSolver.solverStep(time, dt)
            time += dt

            #Collect waveField values
            #hdf5.collect("waveField", time, dt)

            i += 1
        #Export pressure at receivers positions to .segy file
        #exportToSegy(table = pressureAtReceivers.to_numpy(),
        #             shot = shot,
        #             filename = "pressure_Shot"+shot.id,
        #             directory = acquisition.output,
        #             rank = rank)

        #Output waveField values
        #hdf5.output("waveField", time, dt)

        shot.flag = "Done"
        resetWaveField(problem)

        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1

    pygeosx._finalize()


if __name__ == "__main__":
    main()

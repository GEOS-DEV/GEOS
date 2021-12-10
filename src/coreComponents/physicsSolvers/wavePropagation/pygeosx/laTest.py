import pygeosx
import sys
from seismicUtilities.acquisition import EQUISPACEDAcquisition
#from seismicUtilities.segy import exportToSegy
from seismicUtilities.acoustic import updateSourceAndReceivers, resetWaveField, setTimeVariables
#from mpi4py import MPI

#comm = MPI.COMM_WORLD
#rank = comm.Get_rank()


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
                                        number_of_receivers = 100,
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
    acquisition.limitedAperture(500)
    acquisition.calculDt()


    maxTime = 2.0
    outputSeismoTraceInterval = 5
    outputWaveFieldInterval = 100


    result = acousticShot(maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acquisition)



def acousticShot(maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acquisition, rank=0):
    #Loop over the shots
    segyList = []
    ishot=0
    print(acquisition)
    for shot in acquisition.shots:
        sys.argv[2] = shot.xml
        dt = shot.dt
        dtSeismoTrace = dt * outputSeismoTraceInterval
        nsamples = int(maxTime/(dt*outputSeismoTraceInterval))+1

        #Initialize GEOSX, get Acoustic Solver, get HDF5 Outputs
        if ishot == 0:
            problem = pygeosx.initialize(rank, sys.argv)
        else:
            problem = pygeosx.reinit(sys.argv)

        acousticSolver = pygeosx.pysolver.Solver("/Solvers/acousticSolver")
        hdf5 = pygeosx.pyhdf5.HDF5()

        #Get view on pressure at receivers locations
        pressureAtReceivers = acousticSolver.get_wrapper("pressureNp1AtReceivers").value()

        #Set time variables and update source and receivers positions
        setTimeVariables(problem, maxTime, dtSeismoTrace)
        updateSourceAndReceivers(acousticSolver, shot.source, shot.receivers)
        pygeosx.apply_initial_conditions()

        shot.flag = "In Progress"

        time = 0.0
        i = 0
        #Time Loop
        while time < maxTime:
            if rank == 0:
                print("time = %.3fs," % time, "dt = %.4f," % dt, "iter =", i+1)

            acousticSolver.explicitStep(time, dt)

            time += dt
            i += 1
            #Collect waveField values
            if i % outputWaveFieldInterval == 0:
                hdf5.collect("waveField", time, dt)

        #exportToSegy(table = pressureAtReceivers.to_numpy(),
        #             shot = shot,
        #             filename = "pressure_Shot"+shot.id,
        #             directory = acquisition.output,
        #             rank = rank)
        segyList.append("pressure_Shot"+shot.id)
        print(pressureAtReceivers.to_numpy())
        #Output waveField values
        hdf5.output("waveField", time, dt)

        shot.flag = "Done"
        resetWaveField(problem)

        if rank == 0:
            print("Shot", shot.id, "done\n")

        ishot += 1

    pygeosx._finalize()

    return segyList


if __name__ == "__main__":
    main()

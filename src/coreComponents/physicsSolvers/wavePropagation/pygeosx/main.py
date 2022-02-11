'''
Created on 9/02/2021

@author: macpro
'''
from seismicUtilities.acquisition import EQUISPACEDAcquisition
from seismicUtilities.client import Client, LSFCluster,SLURMCluster, Future
#from seismicUtilities.AcousticShot import acoustic_shots
from laTest import acousticShot
from utils import parse_args
import time

def main():
    sysargs = parse_args()

    acq = EQUISPACEDAcquisition(boundary=[[0,13520],[0,13520],[0,4200]],
                                dt=0.002,
                                velocity_model="/beegfs/jbesset/codes/SEP_REDUCE_Model/338x338x105_velModel.geos",
                                start_source_pos    = [4001, 6001],
                                end_source_pos      = [11001, 6001],
                                start_receivers_pos = [[51, 6001]],
                                end_receivers_pos   = [[13491, 6001]],
                                number_of_sources   = 2,
                                number_of_receivers = 500,
                                source_depth = 4099,
                                receivers_depth = 4149)


    acq.add_xml(sysargs.xml)
    acq.limitedAperture(2000)
    #acq.calculDt()

    acqs = acq.split(2)

    cluster = SLURMCluster(job_name="seismicAcquisition", nodes=1, cores=4,
                           python="/beegfs/jbesset/codes/GEOSX/build-test_module-release/lib/PYGEOSX/bin/python",
                           launch="mpirun",
                           extra=["-C bora"],
                           env_extra=["physics/geosx_deps/mkl2020_gcc9.3.0_mpi4.0.3/NoCuda"])

    client = Client(cluster)
    client.scale(2)


    maxTime = 2.0
    outputSeismoTraceInterval = 5
    outputWaveFieldInterval = 100


    args = []
    for acq in acqs:
        args.append((maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acq,))

    #args = (maxTime, outputSeismoTraceInterval, outputWaveFieldInterval, acq,)


    futures = client.map(acousticShot, args,
                         cmd_line = ["-i", sysargs.xml],
                         cores=4, x_partition=4)

    #future = client.submit(acousticShot, args,
    #                       cmd_line=["-i", sysargs.xml],
    #                       cores=2, x_partition=2)
if __name__ == "__main__":
    main()

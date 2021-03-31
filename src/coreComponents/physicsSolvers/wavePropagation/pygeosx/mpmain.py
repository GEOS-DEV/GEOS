'''
Created on 9/02/2021

@author: macpro
'''

import pygeosx
import sys
from mpi4py import MPI
import multiprocessing as mp

from mesh import *
from acquisition import *
from shotFileManager import *


def multiProcessing(shot_file):
    cmd = "python /home/m3d/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/pygeosx/main.py -i /home/m3d/codes/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/benchmarks/pygeosx_test.xml " + str(shot_file)

    os.system(cmd)



def main():

    problem = pygeosx.initialize(0, sys.argv)
    pygeosx.apply_initial_conditions()
    maxT = problem.get_wrapper("Events/maxTime").value()[0]
    dt = calculDt(problem)
    boundary_box  = domainBoundary(problem)
    # - This is a dummy initialization to be able to get the seismic acquisition/shot_list and split it between all proc
    # - Note that the current way to get a seismic acquisition will be replaced by the lecture of
    #   a segy file ==> no longer need to do this dummy pygeosx.initialize()


    frequency = 10.0
    wavelet   = ricker(maxT, dt, frequency)

    folder_path = '/home/m3d/Desktop/pygeosx/sismoTrace/'
    shot_list = segy_acquisition(folder_path, wavelet, dt)

    """
    #Set acquisition parameters
    nb_source_x     = 5
    nb_source_y     = 2
    nb_receiver_x   = 3
    nb_receiver_y   = 4
    receiver_zone_x = 50
    receiver_zone_y = 90
    #Get seismic acquisition
    shot_list = moving_acquisition(boundary_box, wavelet, nb_source_x, nb_source_y, nb_receiver_x, nb_receiver_y, receiver_zone_x, receiver_zone_y)
    """

    """At this point we have defined our seismic acquisition/shot_list"""

    nb_proc = 1
    p = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    #Loop over the process launch
    for i in range(nb_proc):
        nb_shot = int(nb_shot_m1/(nb_proc-i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])

        p.append( mp.Process(target = multiProcessing,
                             args   = (shot_file,) ) )

        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

        #Start process
        p[i].start()


    for i in range(nb_proc):
        p[i].join()


if __name__ == "__main__":
    main()

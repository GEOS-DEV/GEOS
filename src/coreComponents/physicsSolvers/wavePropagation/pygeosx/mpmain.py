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


def multiProcessing(shot_file, basePath, xmlPath):

    wavePropagationPath = basePath + "GEOSX/src/coreComponents/physicsSolvers/wavePropagation/"
    pygeosxPath = wavePropagationPath + "/pygeosx/"
    tracePath = pygeosxPath + "/sismoTrace/"

    cmd = "python " + pygeosxPath + "/main.py -i " + xmlPath + " " + shot_file + " " + tracePath

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

    basePath = "/home/m3d/codes/"

    wavePropagationPath = basePath + "/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/"
    pygeosxPath = wavePropagationPath + "/pygeosx/"
    segyPath = pygeosxPath + "/sismoTrace/"
    xmlPath = sys.argv[2]

    shot_list = segy_acquisition(segyPath, wavelet, dt)

    #Get seismic acquisition
    """
    shot_list = moving_acquisition(boundary_box,
                                   wavelet,
                                   dt,
                                   nbsourcesx = 5,
                                   nbsourcesy = 2,
                                   nbreceiversx = 3,
                                   nbreceiversy = 4,
                                   lenRx = 50,
                                   lenRy = 90)
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
                             args   = (shot_file, basePath, xmlPath) ) )

        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

        #Start process
        p[i].start()


    for i in range(nb_proc):
        p[i].join()


if __name__ == "__main__":
    main()

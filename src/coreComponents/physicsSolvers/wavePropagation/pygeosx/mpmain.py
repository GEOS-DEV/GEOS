'''
Created on 9/02/2021

@author: macpro
'''

import pygeosx
import sys
import multiprocessing as mp

from mesh import *
from acquisition import *
from shotFileManager import *

basePath = "/home/m3d/codes/"

wavePropagationPath = basePath + "/GEOSX/src/coreComponents/physicsSolvers/wavePropagation/"
pygeosxPath = wavePropagationPath + "/pygeosx/"
segyPath = pygeosxPath + "/segyAcquisition/"
xmlPath = sys.argv[2]


def multiProcessing(shot_list, nb_proc = 1):
    p = []
    nb_shot_m1 = len(shot_list)
    ind = 0

    #Loop over the process launch
    for i in range(nb_proc):
        nb_shot = int(nb_shot_m1/(nb_proc-i))

        shot_file = exportShotList(i, shot_list[ind:ind + nb_shot])

        p.append( mp.Process(target = mainProcess,
                             args   = (shot_file,) ) )

        ind = ind + nb_shot
        nb_shot_m1 = nb_shot_m1 - nb_shot

        #Start process
        p[i].start()


    for i in range(nb_proc):
        p[i].join()


def mainProcess(shot_file):

    tracePath = os.path.abspath(os.getcwd()) + "/outputSismoTrace/"
    if os.path.exists(tracePath):
        pass
    else:
        os.mkdir(tracePath)

    tracePath = tracePath + "/traceProc" + str(mp.current_process()._identity)[1] + "/"
    if os.path.exists(tracePath):
        pass
    else:
        os.mkdir(tracePath)


    cmd = "mpirun -np 2 python " + pygeosxPath + "/main.py -i " + xmlPath + " -x 2 " + shot_file + " " + tracePath

    os.system(cmd)



def main():

    os.system("python firstInit.py " + str(sys.argv[1]) + " " + str(sys.argv[2]))

    maxT, dt, boundary_box = readInitVariable("init_variable.txt")

    frequency = 10.0
    wavelet   = ricker(maxT, dt, frequency)


    shot_list = equispaced_acquisition(boundary_box,
                                       wavelet,
                                       dt,
                                       [1001, 11001, 4],
                                       [6751, 6751, 1],
                                       13481,
                                       [21, 13481, 5],
                                       [6751, 6751, 1],
                                       13491,
                                       export = 1
                                       )


    """
    shot_list = moving_acquisition(boundary_box,
                                   wavelet,
                                   dt,
                                   nbsourcesx = 30,
                                   nbsourcesy = 1,
                                   nbreceiversx = 50,
                                   nbreceiversy = 1,
                                   lenRx = 50,
                                   lenRy = 0,
                                   export = 1)

    """

    #shot_list = segy_acquisition(segyPath + acqName, wavelet, dt)

    #multiProcessing(shot_list)



if __name__ == "__main__":
    main()

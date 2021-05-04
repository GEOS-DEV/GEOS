'''
Created on 9/02/2021

@author: macpro
'''

import pygeosx
import sys
from processSpawn import firstInit, multiProcessing, daskProcessing
from seismicAcquisition import equispaced_acquisition
from source import ricker

def main():

    maxT, dt = firstInit(sys.argv[3], sys.argv[2], mpiranks = 2, x = 2)

    frequency = 5.0
    wavelet   = ricker(maxT, dt, frequency)


    shot_list = equispaced_acquisition([[0,2000],[0,2000],[0,2000]],
                                       wavelet,
                                       dt,
                                       [101, 1001],
                                       [1901, 1001],
                                       [[21, 1001],[1001, 21]],
                                       [[1981, 1001],[1001, 1981]],
                                       number_of_sources = 1,
                                       number_of_receivers = [10]
                                       )


    #multiProcessing(shot_list, sys.argv[3], sys.argv[2])
    daskProcessing(shot_list, sys.argv[3], sys.argv[2])

if __name__ == "__main__":
    main()

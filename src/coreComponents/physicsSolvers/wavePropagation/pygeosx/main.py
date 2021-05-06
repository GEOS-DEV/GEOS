'''
Created on 9/02/2021

@author: macpro
'''

import pygeosx
from processSpawn import firstInit, multiProcessing, daskProcessing
from seismicAcquisition import equispaced_acquisition
from source import ricker
from fileManager import parse_args

args = parse_args()

def main():

    maxT, dt = firstInit(args.geosx, args.xml, mpiranks = 2, x = 2)

    frequency = 5.0
    wavelet   = ricker(maxT, dt, frequency)


    shot_list = equispaced_acquisition([[0,2000],[0,2000],[0,2000]],
                                       wavelet,
                                       dt,
                                       start_source_pos    = [101, 1001],
                                       end_source_pos      = [1901, 1001],
                                       start_receivers_pos = [[21, 1001], [1001, 21]],
                                       end_receivers_pos   = [[1981, 1001], [1001, 1981]],
                                       number_of_sources   = 2,
                                       number_of_receivers = [10]
                                       )

    #multiProcessing(shot_list, args.geosx, args.xml)
    #daskProcessing(shot_list, args.geosx, args.xml)

if __name__ == "__main__":
    main()

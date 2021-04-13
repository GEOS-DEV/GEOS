import pygeosx
from mpi4py import MPI
import sys

from shotSimulation import *
from shotFileManager import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():

    shot_file = sys.argv[5]
    tracePath = sys.argv[6]

    shot_list = readShotList(shot_file)

    shot_simul(rank, sys.argv[0:5], shot_list, tracePath)


if __name__ == "__main__":
    main()

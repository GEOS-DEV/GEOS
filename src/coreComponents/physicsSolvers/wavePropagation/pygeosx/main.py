import pygeosx
from mpi4py import MPI
import sys

from shotSimulation import *
from shotFileManager import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():

    shot_file = sys.argv[7]
    
    tracePath = sys.argv[8]

    shot_list = readShotList(shot_file)

    shot_simul(rank, sys.argv[0:7], shot_list, tracePath)

    if rank==0:
        path = os.path.abspath(os.getcwd()) + "/shots_lists/"
        os.remove(shot_file)
        os.rmdir(path)

if __name__ == "__main__":
    main()

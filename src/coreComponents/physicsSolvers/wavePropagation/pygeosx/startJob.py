from mpi4py import MPI
import sys
import pygeosx
from utils import *
import importlib

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

acquisition = json_to_obj(sys.argv[3])
module = importlib.import_module(sys.argv[1])
func = eval(sys.argv[2])

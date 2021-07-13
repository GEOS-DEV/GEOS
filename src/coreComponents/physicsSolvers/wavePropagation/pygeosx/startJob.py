from mpi4py import MPI
import sys
import pygeosx
from utils import *
import time
import importlib

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

json = open(sys.argv[11])
acquisition = json_to_obj(json)
json.close()

if rank==0:
    time.sleep(1)
    os.remove(sys.argv[11])

module_str = sys.argv[9]
func_str = sys.argv[10]

module = importlib.import_module(module_str)
func = getattr(module, func_str)

problem = pygeosx.initialize(rank, sys.argv[0:9])
pygeosx.apply_initial_conditions()
func(rank, problem, acquisition, acquisition.output)

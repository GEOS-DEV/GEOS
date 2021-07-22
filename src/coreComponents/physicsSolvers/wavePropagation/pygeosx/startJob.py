from mpi4py import MPI
import sys
import pygeosx
from utils import *
import time
import importlib
from acquisition import Acquisition

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

json = open(sys.argv[11])
dic = json_to_dict(json)
json.close()
acquisition = Acquisition().construct_from_dict(**dic)

print(acquisition.shot)

if rank==0:
    time.sleep(1)
    os.remove(sys.argv[11])

    outputDir = sys.argv[12]
    keyFile = sys.argv[13]+".txt"
    outputFile = os.path.join(outputDir, keyFile)


module_str = sys.argv[9]
func_str = sys.argv[10]

module = importlib.import_module(module_str)
func = getattr(module, func_str)

problem = pygeosx.initialize(rank, sys.argv[0:9])
pygeosx.apply_initial_conditions()
result = func(rank, problem, acquisition)

if rank == 0:
    with open(outputFile, 'w') as f:
        for line in result:
            f.write(line + "\n")
    f.close()

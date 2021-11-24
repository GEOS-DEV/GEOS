from mpi4py import MPI
import sys
import os
sys.path.append(os.getcwd() + '/..')
import pygeosx
from utils import *
import time
import importlib
from acquisition import Acquisition

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


with open(sys.argv[9], 'r') as f:
    dict_args = json.load(f)
f.close()
args = []

for k, v in dict_args.items():
    if not isinstance(v, dict):
        args.append(v)
    else:
        args.append(munchify(v))

"""
json = open(sys.argv[9])
dict = json_to_dict(json)
json.close()

acquisition = Acquisition().construct_from_dict(**dict)
"""

xml = args[0].shots[0].xml

if rank==0:
    time.sleep(1)
    os.remove(sys.argv[9])

    outputDir = sys.argv[10]
    keyFile = sys.argv[11]+".txt"
    outputFile = os.path.join(outputDir, keyFile)


module_str = sys.argv[7]
func_str = sys.argv[8]

module = importlib.import_module(module_str)
func = getattr(module, func_str)

geosx_argv = [sys.argv[0], "-i", xml]
for arg in sys.argv[1:7]:
    geosx_argv.append(arg)

problem = pygeosx.initialize(rank, geosx_argv)
pygeosx.apply_initial_conditions()
result = func(*args)

if rank == 0:
    with open(outputFile, 'w') as f:
        for line in result:
            f.write(line + "\n")
    f.close()

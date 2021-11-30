from mpi4py import MPI
import sys
import os
import pygeosx
import time
import importlib
import json
from munch import munchify


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


with open(sys.argv[-3], 'r') as f:
    dict_args = json.load(f)
f.close()
args = []


if rank==0:
    time.sleep(1)
    os.remove(sys.argv[-3])

    outputDir = sys.argv[-2]
    keyFile = sys.argv[-1]+".txt"
    outputFile = os.path.join(outputDir, keyFile)

module_str = sys.argv[-5]
func_str = sys.argv[-4]


for k, v in dict_args.items():
    if not isinstance(v, dict):
        args.append(v)
    else:
        args.append(munchify(v))

geosx_argv = [sys.argv[0], "-i", ""]
for arg in sys.argv[1:-5]:
    geosx_argv.append(arg)

args.append(geosx_args)
args.append(rank)

module = importlib.import_module(module_str)
func = getattr(module, func_str)

result = func(*args)

if rank == 0:
    with open(outputFile, 'w') as f:
        for line in result:
            f.write(line + "\n")
    f.close()

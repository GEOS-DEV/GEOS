from mpi4py import MPI
import sys
import os
sys.path.append(os.getcwd())
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

cmd_args = [sys.argv[0]]
for arg in sys.argv[1:-5]:
    cmd_args.append(arg)

sys.argv = cmd_args

module = importlib.import_module(module_str)
func = getattr(module, func_str)

result = func(*args, rank=rank)

if rank == 0:
    with open(outputFile, 'w') as f:
        for line in result:
            f.write(line + "\n")
    f.close()

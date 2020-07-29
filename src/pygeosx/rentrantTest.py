import pygeosx
import sys
from mpi4py import MPI

# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
  print( "In python" )
  sys.stdout.flush()

problem = pygeosx.initialize( rank, sys.argv )
pygeosx.applyInitialConditions()

while pygeosx.run() != pygeosx.COMPLETED:
    pass

if rank == 0:
  print("\n\n In python second time around")

problem = pygeosx.reinit( sys.argv )
pygeosx.applyInitialConditions()

while pygeosx.run() != pygeosx.COMPLETED:
    pass

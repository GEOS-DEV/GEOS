import pygeosx
import sys
from mpi4py import MPI

# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
    print( "In python" )
    sys.stdout.flush()

# If not a restart run then we'll run the SSLE-sedov problem first. 
if "-r" not in sys.argv:
  argCopy = list( sys.argv )

  xmlIndex = argCopy.index( "-i" ) + 1
  argCopy[ xmlIndex ] = "/usr/WS2/corbett5/geosx/integratedTests/update/run/solidMechanicsSSLE/SSLE-sedov.xml"

  outputIndex = argCopy.index( "-o" ) + 1
  argCopy[ outputIndex ] = sys.argv[ outputIndex ] + "_ssle_dummy"

  problem = pygeosx.initialize( rank, argCopy )
  pygeosx.applyInitialConditions()

  while pygeosx.run() != pygeosx.COMPLETED:
      pass

  if rank == 0:
    print("\n\n In python second time around")
    sys.stdout.flush()

  problem = pygeosx.reinit( sys.argv )
else:
  problem = pygeosx.initialize( rank, sys.argv )

pygeosx.applyInitialConditions()

while pygeosx.run() != pygeosx.COMPLETED:
    pass

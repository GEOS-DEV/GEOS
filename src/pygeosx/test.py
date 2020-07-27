import pygeosx
import sys
from mpi4py import MPI

# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def printAndFlush( msg ):
    print( "Rank {}: {}".format( rank, msg ) )
    sys.stdout.flush()



printAndFlush( "In python before initialization." )

pygeosx.initialize( rank, sys.argv )
currentTime = pygeosx.get( "Events/time", False )
printAndFlush( "In python after initialization: current time = {}".format( currentTime[ 0 ] ) )

pygeosx.applyInitialConditions()
currentTime = pygeosx.get( "Events/time", False )
printAndFlush( "In python after applyInitialConditions: current time = {}".format( currentTime[ 0 ] ) )

while pygeosx.run() != pygeosx.COMPLETED:
    currentTime = pygeosx.get( "Events/time", True )
    printAndFlush( "In python: current time = {}".format( currentTime[ 0 ] ) )
    currentTime[ 0 ] += 1e-6

currentTime = pygeosx.get( "Events/time", False )
printAndFlush( "In python after after the simulation has ended: current time = {}".format( currentTime[ 0 ] ) )

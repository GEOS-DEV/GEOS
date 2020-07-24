import pygeosx
import sys
import os

xml = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), "pyssle.xml" )

print( "In python before initialization." )
sys.stdout.flush()

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_ranks = comm.Get_size()

# rank = 0

pygeosx.initialize( rank, sys.argv )
currentTime = pygeosx.get( "Events/time", False )
print( "In python after initialization: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

pygeosx.applyInitialConditions()
currentTime = pygeosx.get( "Events/time", False )
print( "In python after applyInitialConditions: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

while pygeosx.run() != pygeosx.COMPLETED:
    currentTime = pygeosx.get( "Events/time", True )
    print( "In python: current time = {}".format( currentTime[ 0 ] ) )
    sys.stdout.flush()
    currentTime[ 0 ] += 1e-6

currentTime = pygeosx.get( "Events/time", False )
print( "In python after after the simulation has ended: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

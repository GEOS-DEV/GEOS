import pygeosx
import sys
import os
from mpi4py import MPI

xml = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), "pyssle.xml" )

print( "In python before initialization." )
sys.stdout.flush()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
n_ranks = comm.Get_size()

pygeosx.initialize( rank, sys.argv )
currentTime = pygeosx.get( "Events/time", False )
print( "In python after initialization: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

pygeosx.applyInitialConditions()
currentTime = pygeosx.get( "Events/time", False )
print( "In python after applyInitialConditions: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

# 3 is the value of the State::COMPLETED enum
while pygeosx.run() != 3:
    currentTime = pygeosx.get( "Events/time", True )
    print( "In python: current time = {}".format( currentTime[ 0 ] ) )
    sys.stdout.flush()
    currentTime[ 0 ] += 1e-6

currentTime = pygeosx.get( "Events/time", False )
print( "In python after after the simulation has ended: current time = {}".format( currentTime[ 0 ] ) )
sys.stdout.flush()

pygeosx.finalize()

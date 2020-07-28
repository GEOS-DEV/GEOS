import pygeosx
import sys
from mpi4py import MPI

# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def printAndFlush( msg ):
    print( "Rank {}: {}".format( rank, msg ) )
    sys.stdout.flush()

def printWithIndent( msg, indent ):
    indentString = " " * indent
    print( indentString + msg.replace( "\n", "\n" + indentString ) )

def printGroup( group, indent=0 ):
    print( "{}{}".format( " " * indent, group ) )

    indent += 4
    print( "{}wrappers:".format( " " * indent ) )

    for wrapper in group.wrappers():
        print( "{}{}".format( " " * ( indent + 4 ), wrapper ) )
        printWithIndent( str( wrapper.value( False ) ), indent + 8 )

    print( "{}groups:".format( " " * indent ) )

    for subGroup in group.groups():
        printGroup( subGroup, indent + 4 )


printAndFlush( "In python before initialization." )

problem = pygeosx.initialize( rank, sys.argv )

currentTime = problem.getWrapper( "Events/time" ).value( False )
printAndFlush( "In python after initialization: current time = {}".format( currentTime[ 0 ] ) )

pygeosx.applyInitialConditions()
currentTime = problem.getWrapper( "Events/time" ).value( False )
printAndFlush( "In python after applyInitialConditions: current time = {}".format( currentTime[ 0 ] ) )

while pygeosx.run() != pygeosx.COMPLETED:
    currentTime = problem.getWrapper( "Events/time" ).value( True )
    printAndFlush( "In python: current time = {}".format( currentTime[ 0 ] ) )
    currentTime[ 0 ] += 1e-6

currentTime = problem.getWrapper( "Events/time" ).value( False )
printAndFlush( "In python after after the simulation has ended: current time = {}".format( currentTime[ 0 ] ) )

printGroup( problem.getGroup( "domain/MeshBodies/mesh1/Level0/nodeManager" ) )



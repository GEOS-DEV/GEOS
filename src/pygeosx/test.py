import pygeosx
import sys
import os

xml = os.path.join( os.path.dirname( os.path.realpath( __file__ ) ), "pyssle.xml" )

print( "In python before initialization." )

pygeosx.initialize( [sys.argv[ 0 ], "-i", xml] )

currentTime = pygeosx.get( "Events/time", False )
print( "In python: current time = {}".format( currentTime[ 0 ] ) )

while pygeosx.run():
    currentTime = pygeosx.get( "Events/time", True )
    print( "In python: current time = {}".format( currentTime[ 0 ] ) )
    currentTime[ 0 ] += 1e-6

pygeosx.finalize()

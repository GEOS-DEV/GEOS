import os
import sys
import argparse
import re
import subprocess
import json
import glob

class style():
    RED = '\033[31m'
    GREEN = '\033[32m'
    RESET = '\033[0m'

resultRegex = r"init time = (.*)s, run time = (.*)s"
#caliQueryBin = "/usr/gapps/GEOSX/thirdPartyLibs/2020-10-15/install-lassen-clang@upstream-release/caliper/bin/cali-query"
caliQueryBin = "/usr/gapps/GEOSX/thirdPartyLibs/2020-10-15/install-quartz-clang@10.0.0-release/caliper/bin/cali-query"
caliQueryFString = ' {} -q "SELECT * GROUP BY function FORMAT json ORDER BY time.duration" '
#caliQueryString = ' {} -q "SELECT * WHERE function={} GROUP BY function FORMAT json ORDER BY time.duration" '

def getTimesFromFile( filePath ):
    """
    Return the init time and run time from a GEOSX standard output file.

    Arguments:
        filePath: The path of the output file to parse.
    """
    with open( filePath, "r" ) as file:
        for line in file:
            matches = re.search( resultRegex, line )
            if matches is not None:
                return [ float( matches.groups()[ 0 ] ), float( matches.groups()[ 1 ] ) ]


def getTimesFromFolder( folder ):
    """
    Arguments:
        folder: The top level directory the benchmarks were run in.
    """
    results = {}
    for outerFile in os.listdir( folder ):
        xmlName = outerFile
        outerFile = os.path.join( folder, outerFile )

        if os.path.isdir( outerFile ):
            for innerFile in os.listdir( outerFile ):
                problemName = innerFile
                innerFile = os.path.join( outerFile, innerFile )

                if os.path.isdir( innerFile ):
                    outputFile = os.path.join( innerFile, "output.txt" );

                    if not os.path.exists( outputFile ) or not os.path.isfile( outputFile ):
                        print( "{} does not exists or is not a file!".format( outputFile ) )

                    init_run = getTimesFromFile( outputFile )
                    if init_run is not None:
                        nodeCount = parseNodesFromProblemName( problemName )
                        # print( nodeCount )

                        if results.get(nodeCount, None) is None:
                            results[nodeCount] = {}
                        elemCount = parseElementsFromFilename( xmlName )
                        # print( elemCount )

                        if results[nodeCount].get(elemCount, None) is None:
                            results[nodeCount][elemCount] = { }
                        results[nodeCount][elemCount]["init"] = init_run[0]
                        results[nodeCount][elemCount]["run"] = init_run[1]

                    caliFile = glob.glob( os.path.join( innerFile, "*.cali" ) )
                    if len(caliFile) > 0:
                        caliFile = caliFile[0]
                        if os.path.exists( caliFile ):
                            caliQuery = caliQueryBin + caliQueryFString.format( caliFile )
                            #print( caliQuery )
                            caliJson = subprocess.check_output( caliQuery, shell=True )
                            caliData = json.loads( caliJson )
                            for entry in caliData:
                                funcPath = entry.get("path", "")
                                time = entry.get("max#inclusive#sum#time.duration", "")
                                results[nodeCount][elemCount][funcPath] = float( time )
    return results

elemRe = re.compile( ".*\_([0-9]+)M_elem" )
def parseElementsFromFilename( xmlName ):
    # print( xmlName )
    matches = elemRe.match( xmlName )
    millionsElements = matches.groups()[ 0 ]
    return int( millionsElements )

nodesRe = re.compile( "(?:MPI_)*(?:OMP_)*(?:CUDA_)*([0-9]*)" )
def parseNodesFromProblemName( problemName ):
    # print( problemName )
    matches = nodesRe.match( problemName )
    nodes = matches.groups()[ 0 ]
    return int( nodes )

def generateTable( results, measure ):
    """
    Print a table containing the speed up of the results over the baseline results.

    Arguments:
        results: The dictionary of scaling results.
    """

    nodeCounts = results.keys()
    elemCounts = []
    for nodeCount in nodeCounts:
        elemCounts.extend( results[nodeCount].keys() )
    elemCounts = list( set( elemCounts ) )

    strings = [ "{:15}".format("ElemCol/NodeRow") ]
    for elemCount in sorted(elemCounts):
        strings += [ "{:15}".format( elemCount ) ]
    strings += [ "\n" ]

    for nodeCount in sorted(nodeCounts):
        strings += [ "{:15}".format( nodeCount ) ]
        for elemCount in sorted(elemCounts):
            result = results[nodeCount].get(elemCount)
            if result is not None:
                val = result.get( measure, " " )
                strings += [ "{:15}".format( val ) ]
        strings += [ "\n" ]

    print( ",".join( strings ) )

def printMeasures( results ):
    nodeCounts = results.keys()
    elemCounts = []
    for nodeCount in nodeCounts:
        elemCounts.extend( results[nodeCount].keys() )
    elemCounts = list( set( elemCounts ) )

    strings = []
    for nodeCount in sorted(nodeCounts):
        for elemCount in sorted(elemCounts):
            result = results[nodeCount].get(elemCount)
            if result is not None:
                strings.extend( result.keys( ) )

    strings = sorted( list( set( strings ) ) )
    print( "\n".join( strings ) )


def main():
    """ Parse the command line arguments and compare the benchmarks. """

    parser = argparse.ArgumentParser()
    parser.add_argument( "scalingDir", help="The directory where the scaling study was run" )
    parser.add_argument( "--measure", help="What measure to print scaling data for (init/run)" )
    parser.add_argument( "--list-measures", action='store_const', const=True, help="List all the available measures.")
    args = parser.parse_args()

    scalingDir = os.path.abspath( args.scalingDir )
    if not os.path.isdir( scalingDir ):
        raise ValueError( "scalingDir is not a directory!" )

    results = getTimesFromFolder( scalingDir )

    if args.list_measures:
        printMeasures( results )

    if args.measure:
        generateTable( results, args.measure )

    return 0


if __name__ == "__main__" and not sys.flags.interactive:
    sys.exit(main())

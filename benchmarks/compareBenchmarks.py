import os
import sys
import argparse
import re


class style():
    RED = '\033[31m'
    GREEN = '\033[32m'
    RESET = '\033[0m'


resultRegex = r"init time = (.*)s, run time = (.*)s"


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
                return float( matches.groups()[ 0 ] ), float( matches.groups()[ 1 ] )
    
    raise Exception( "Could not get times from {}".format( filePath ) )


def getTimesFromFolder( folder ):
    """
    Return a dictionary containing the init and run times of each run in the benchmark folder.

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
                        raise ValueError( "{} does not exist or is not a file.".format( outputFile ) )

                    init, run = getTimesFromFile( outputFile )
                    results[ ( xmlName, problemName ) ] = init, run
    
    return results


def joinResults( results, baselineResults ):
    """
    Return a dictionary containing both the results and baseline results.
    
    Arguments:
        results: The dictionary of benchmark results.
        baselineResults: The dictionary of baseline benchmark results.
    """
    joined = {}
    for key in results:
        joined[ key ] = [ results[ key ][ 0 ], results[ key ][ 1 ], float( "nan" ), float( "nan" ) ]
    
    for key in baselineResults:
        if key in joined:
            joined[ key ][ 2 ] = baselineResults[ key ][ 0 ]
            joined[ key ][ 3 ] = baselineResults[ key ][ 1 ]
        else:
            joined[ key ] = [ float( "nan" ), float( "nan" ), baselineResults[ key ][ 0 ], baselineResults[ key ][ 1 ] ]
    
    joinedList = []
    for key in joined:
        item = []
        item += key
        item += joined[ key ]
        joinedList.append( item )
    
    return sorted( joinedList, lambda x, y: cmp( x, y ) )


def getValue( x ):
    """
    If x is a tuple return the first entry, else return x.

    Arguments:
        x: The object to get the value of.
    """
    if isinstance( x, tuple ):
        return x[ 0 ]
    else:
        return x


def getColor( x ):
    """
    If x is a tuple return the second entry, which should be an ANSI color code. Else return the default color.

    Arguments:
        x: The object to get the color of.
    """
    if isinstance( x, tuple ):
        return x[ 1 ]
    else:
        return style.RESET



def printTable( table ):
    """
    Print a table in a nice format, with optional coloring.

    Arguments:
        table: A list of rows to print. Each row should be of the same length. Then entries in each row
            should either be a string or a tuple of a string and ANSI color code.
    """
    col_width = [ max( len( getValue( x ) ) for x in col ) for col in zip( *table ) ]
    print( "| " + " | ".join( "{:{}}".format( getValue( x ), col_width[ i ] ) for i, x in enumerate( table[ 0 ] ) ) + " |" )
    print( "|" + "|".join( "-" * width + "--" for width in col_width ) + "|" )

    for line in table[ 1: ]:
        print( "| " + " | ".join( "{}{:{}}{}".format( getColor( x ), getValue( x ), col_width[ i ], style.RESET ) for i, x in enumerate( line ) ) + " |" )

    print( "|" + "|".join( "-" * width + "--" for width in col_width ) + "|" )


def generateTable( results, baselineResults ):
    """
    Print a table containing the speed up of the results over the baseline results.
    
    Arguments:
        results: The dictionary of benchmark results.
        baselineResults: The dictionary of baseline benchmark results.
    """
    lines = [ ( "XML Name", "Problem Name", "init speed up", "run speed up" ) ]
    
    joined = joinResults( results, baselineResults )
    for result in joined:
        xmlName = result[ 0 ]
        problemName = result[ 1 ]
        initTime = result[ 2 ]
        runTime = result[ 3 ]
        baseInitTime = result[ 4 ]
        baseRunTime = result[ 5 ]

        lines.append( ( xmlName, problemName,
                        "{:.2f}x".format( baseInitTime / initTime ),
                        "{:.2f}x".format( baseRunTime / runTime ) ) )

    printTable( lines )


def main():
    """ Parse the command line arguments and compare the benchmarks. """

    parser = argparse.ArgumentParser()
    parser.add_argument( "toCompareDir", help="The directory where the new benchmarks were run." )
    parser.add_argument( "baselineDir", help="The directory where the baseline benchmarks were run." )
    args = parser.parse_args()

    toCompareDir = os.path.abspath( args.toCompareDir )
    if not os.path.isdir( toCompareDir ):
        raise ValueError( "toCompareDir is not a directory!" )
    
    baselineDir = os.path.abspath( args.baselineDir )
    if not os.path.isdir( baselineDir ):
        raise ValueError( "baselineDir is not a directory!" )

    results = getTimesFromFolder( toCompareDir )
    baselineResults = getTimesFromFolder( baselineDir )

    generateTable( results, baselineResults )
    return 0


if __name__ == "__main__" and not sys.flags.interactive:
    sys.exit(main())

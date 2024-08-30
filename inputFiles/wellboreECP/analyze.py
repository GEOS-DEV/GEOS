import os
import sys
import argparse
import re

class style():
    RED = '\033[31m'
    GREEN = '\033[32m'
    RESET = '\033[0m'

def getMetricsFromFile( filePath, **regexes ):
    """
    Return matching metrics from an output file.

    Arguments:
        filePath: The path of the output file to parse.
    """
    with open( filePath, "r" ) as file:
        data = file.read() 
    
    times = {}
    for title, regex in regexes.items():
        matches = regex.search( data )
        if matches is not None:
            metric = tuple( matches.groups() )
            if len(metric) == 1:
                metric = metric[0]
            times[title] = metric
    if len(times) == 0:
        raise ValueError( "Could not get times from {}".format( filePath ) )
    return times

def getMetricsFromDirectory( problemDir, filere, **regexes ):
    """
    Return a dictionary containing the regex-described metrics for matching output files in the results directory.

    Arguments:
        problemDir: The top level directory the benchmarks were run in.
    """
    results = {}
    problemName = os.path.basename( problemDir )

    instances = list( filter( lambda filename: os.path.isdir(os.path.join(problemDir, filename)), os.listdir(problemDir) ) )
    for problemInstance in instances:
        instancePath = os.path.join( problemDir, problemInstance )
        matching_files = list( filter( lambda filename: os.path.isfile(os.path.join(instancePath, filename)) and filere.match(filename), os.listdir(instancePath) ) )
        for outputFile in matching_files:
            absOutputFile = os.path.join(instancePath, outputFile)
            try:
                times = getMetricsFromFile( absOutputFile, **regexes )
                results[ ( problemName, problemInstance, outputFile ) ] = times
            except ValueError:
                continue

    return results

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

def printShellTable( table ):
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

def printCsvTable( table, report_dropped ):
    for line in table:
        if len(line) == len(table[0]):
            print(", ".join( str(val) for val in line ) )
    if report_dropped:
        print("\nDropped:")
        print( ", ".join( str(val) for val in table[0] ) )
        for line in table:
            if len(line) != len(table[0]):
                print(", ".join( str(val) for val in line ) )


def generateTable( results, **metrics ):
    """
    Print a table containing the speed up of the results over the baseline results.
    
    Arguments:
        results: The dictionary of benchmark results.
    """
    lines = [ ( "Problem Name", "Instance", "File", *list( metrics.keys( ) ) ) ]
    
    for problem_instance, results in results.items():
        problem, instance, filename = problem_instance
        lines.append( ( problem, instance, filename, *list( results.values( ) ) ) )

    return lines


def analyze( results_dir, filere, metrics, report_dropped ):
    results = getMetricsFromDirectory( results_dir, filere, **metrics )
    table_lines = generateTable( results, **metrics )
    printCsvTable( table_lines, report_dropped )
    return 0

def main( args ):
    parser = argparse.ArgumentParser()
    parser.add_argument( '--results-dir', help="The directory where the benchmarks were run." )
    parser.add_argument( '--filere', help="A filename regex to use to identify results files to parse." )
    parser.add_argument( '--metric', nargs=2, action='append', metavar=('title', 'regex'), help='a parsing metric described ass a title and regex' )
    parser.add_argument( '--report-dropped', action='store_true', default = False, help='print partially collected lines after the full table prints')
    args = parser.parse_args( args )

    results_dir = os.path.abspath( args.results_dir )
    if not os.path.isdir( results_dir ):
        print( f"Error: {results_dir} is not a directory!" )
        sys.exit(1)

    try:
        filere = re.compile( args.filere )
    except re.error:
        print(f"Error: invalid regex '{args.filere}' for file specification.")
        sys.exit(1)

    metrics = {}
    for title, regex_str in args.metric:
        try:
            if regex_str.startswith('"') and regex_str.endswith('"'):
                regex_str = regex_str[1:-1]
            regex = re.compile(regex_str)
            metrics[title] = regex
        except re.error:
            print(f"Error: invalid regex '{regex_str}' for metric '{title}'")
            sys.exit(1)

    return analyze( results_dir, filere, metrics, args.report_dropped )

if __name__ == "__main__" and not sys.flags.interactive:
    args = sys.argv[1:]
    # args = [ "--results-dir", "mechanics", "--filere", "frontier-.*", "--metric", "Git Hash", "sha1:\s*(.*)\)","--metric", "Date", "(.*)\nDone", "--metric", "Run Time", "run time\s*(.*)s", "--metric", "Init Time", "initialization time\s*(.*)s" ]
    sys.exit( main(args) )

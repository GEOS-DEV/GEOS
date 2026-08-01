#!/usr/bin/env python3
"""Harvest and compare linear-solver iteration counts from GEOS log files.

Used to assert iteration-count parity between two runs of the same test suite,
e.g. the hypredrive and legacy hypre solver paths of a hypredrive-enabled build.
Restart checks compare solution fields only; this tool catches solver-quality
regressions (e.g. a preconditioner change that doubles the iteration count while
still converging to the same answer).

Requires the tests to run with ``logLevel >= 1`` on their ``LinearSolverParameters``
so that GEOS emits per-solve lines of the form::

    Linear Solver | Success | Unknowns: ... | Iterations: 12 | ...

Typical usage::

    # After a test-suite run, collect counts keyed by log file path
    compareLinearSolverIterations.py harvest integratedTests/workingDir -o hypredrive.json

    # After re-running the suite through the other solver path
    compareLinearSolverIterations.py harvest integratedTests/workingDir -o legacy.json

    # Compare the two (nonzero exit code on violation)
    compareLinearSolverIterations.py compare hypredrive.json legacy.json
"""

import argparse
import json
import os
import re
import sys

SOLVE_PATTERN = re.compile( r"Linear Solver \| (\w+) \|.*?\| Iterations: (\d+) \|" )

# Log files worth scanning: ATS step logs and plain captures
LOG_SUFFIXES = ( ".log", ".out", ".txt", ".data" )


def harvestFile( path ):
    """Return the per-solve iteration sequence found in a log file (may be empty)."""
    iterations = []
    statuses = []
    try:
        with open( path, errors="replace" ) as f:
            for line in f:
                m = SOLVE_PATTERN.search( line )
                if m:
                    statuses.append( m.group( 1 ) )
                    iterations.append( int( m.group( 2 ) ) )
    except OSError:
        return None
    if not iterations:
        return None
    return { "iterations": iterations,
             "numSolves": len( iterations ),
             "totalIterations": sum( iterations ),
             "maxIterations": max( iterations ),
             "numFailedSolves": sum( 1 for s in statuses if s != "Success" ) }


def harvest( args ):
    results = {}
    root = os.path.abspath( args.directory )
    for dirpath, _, filenames in os.walk( root ):
        for filename in filenames:
            if not filename.endswith( LOG_SUFFIXES ):
                continue
            path = os.path.join( dirpath, filename )
            data = harvestFile( path )
            if data:
                key = os.path.relpath( path, root )
                results[ key ] = data
    with open( args.output, "w" ) as f:
        json.dump( results, f, indent=2, sort_keys=True )
    print( f"Harvested iteration counts from {len( results )} log file(s) under {root}" )
    if not results:
        print( "Warning: no per-solve iteration lines found; "
               "are the tests running with LinearSolverParameters logLevel >= 1?" )
    return 0


def compare( args ):
    with open( args.first ) as f:
        first = json.load( f )
    with open( args.second ) as f:
        second = json.load( f )

    common = sorted( set( first ) & set( second ) )
    onlyFirst = sorted( set( first ) - set( second ) )
    onlySecond = sorted( set( second ) - set( first ) )

    if not common:
        print( "Error: no common log files between the two harvests." )
        return 1

    violations = []
    header = f"{'log file':60s} {'solves':>13s} {'total its':>15s} {'max its':>13s}"
    print( header )
    print( "-" * len( header ) )
    for key in common:
        a, b = first[ key ], second[ key ]
        marks = []

        totalA, totalB = a[ "totalIterations" ], b[ "totalIterations" ]
        relDiff = abs( totalA - totalB ) / max( 1.0, float( max( totalA, totalB ) ) )
        if relDiff > args.total_rel_tol and abs( totalA - totalB ) > args.total_abs_tol:
            marks.append( f"total iterations differ by {100.0 * relDiff:.0f}%" )

        if abs( a[ "maxIterations" ] - b[ "maxIterations" ] ) > args.max_solve_abs_tol:
            marks.append( "worst-solve iteration counts differ by more than "
                          f"{args.max_solve_abs_tol}" )

        if a[ "numFailedSolves" ] != b[ "numFailedSolves" ]:
            marks.append( f"failed-solve counts differ "
                          f"({a[ 'numFailedSolves' ]} vs {b[ 'numFailedSolves' ]})" )

        status = "FAIL: " + "; ".join( marks ) if marks else "ok"
        print( f"{key:60s} {a[ 'numSolves' ]:>6d}/{b[ 'numSolves' ]:<6d} "
               f"{totalA:>7d}/{totalB:<7d} "
               f"{a[ 'maxIterations' ]:>6d}/{b[ 'maxIterations' ]:<6d}  {status}" )
        if marks:
            violations.append( key )

    if onlyFirst or onlySecond:
        print( f"\nNote: {len( onlyFirst )} log file(s) only in {args.first}, "
               f"{len( onlySecond )} only in {args.second} (not compared)." )

    print( f"\nCompared {len( common )} log file(s): "
           f"{len( common ) - len( violations )} ok, {len( violations )} violation(s)." )
    return 1 if violations else 0


def main():
    parser = argparse.ArgumentParser( description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter )
    sub = parser.add_subparsers( dest="command", required=True )

    h = sub.add_parser( "harvest", help="collect iteration counts from log files" )
    h.add_argument( "directory", help="directory to scan recursively" )
    h.add_argument( "-o", "--output", required=True, help="output JSON file" )
    h.set_defaults( func=harvest )

    c = sub.add_parser( "compare", help="compare two harvested JSON files" )
    c.add_argument( "first" )
    c.add_argument( "second" )
    c.add_argument( "--total-rel-tol", type=float, default=0.1,
                    help="allowed relative difference in total iterations per log (default 0.1)" )
    c.add_argument( "--total-abs-tol", type=int, default=5,
                    help="absolute total-iteration slack below which the relative "
                         "criterion is not applied (default 5)" )
    c.add_argument( "--max-solve-abs-tol", type=int, default=5,
                    help="allowed difference in the worst per-solve iteration count (default 5)" )
    c.set_defaults( func=compare )

    args = parser.parse_args()
    return args.func( args )


if __name__ == "__main__":
    sys.exit( main() )

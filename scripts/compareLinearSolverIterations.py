#!/usr/bin/env python3
"""Harvest and compare linear-solver iteration counts from GEOS log files.

Used to assert iteration-count parity between two runs of the same test suite,
e.g. the hypredrive and legacy hypre solver paths of a hypredrive-enabled build.
Restart checks compare solution fields only; this tool catches solver-quality
regressions (e.g. a preconditioner change that shifts a single Newton solve by
one iteration while still converging to the same answer).

Requires the tests to run with ``logLevel >= 1`` on their ``LinearSolverParameters``
so that GEOS emits per-solve lines of the form::

    Linear Solver | Success | Unknowns: ... | Iterations: 12 | ...

Typical usage for the ``_iterative`` hypredrive equivalence family::

    compareLinearSolverIterations.py harvest integratedTests/TestResults/test_data \\
        --iterative -o hypredrive.json
    compareLinearSolverIterations.py harvest integratedTests/TestResults/test_data \\
        --iterative -o legacy.json
    compareLinearSolverIterations.py compare hypredrive.json legacy.json --exact-sequence
"""

import argparse
import json
import os
import re
import sys
import tempfile
import unittest

SOLVE_PATTERN = re.compile( r"Linear Solver \| (\w+) \|.*?\| Iterations: (\d+) \|" )

# ATS prefixes each retained log with the job number, e.g. 0135.geosx_foo_01_1_geos_.log
ATS_JOB_PREFIX = re.compile( r"^\d+\." )

# Log files worth scanning: ATS step logs and plain captures
LOG_SUFFIXES = ( ".log", ".out", ".txt", ".data" )

ITERATIVE_TOKEN = "_iterative"


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


def isGeosxLog( filename ):
    lower = filename.lower()
    return "geosx" in lower and "restartcheck" not in lower


def canonicalKey( relpath, stripAtsPrefix ):
    relpath = relpath.replace( "\\", "/" )
    dirname, basename = os.path.split( relpath )
    if stripAtsPrefix:
        basename = ATS_JOB_PREFIX.sub( "", basename )
    if dirname:
        return f"{dirname}/{basename}"
    return basename


def shouldHarvest( relpath, filename, args ):
    if args.iterative and ITERATIVE_TOKEN not in relpath.replace( "\\", "/" ):
        return False
    if args.geosx_only and not isGeosxLog( filename ):
        return False
    return True


def harvest( args ):
    results = {}
    root = os.path.abspath( args.directory )
    for dirpath, _, filenames in os.walk( root ):
        for filename in filenames:
            if not filename.endswith( LOG_SUFFIXES ):
                continue
            path = os.path.join( dirpath, filename )
            relpath = os.path.relpath( path, root )
            if not shouldHarvest( relpath, filename, args ):
                continue
            data = harvestFile( path )
            if data:
                key = canonicalKey( relpath, args.strip_ats_prefix )
                if key in results:
                    print( f"Error: harvest key collision for {key} "
                           f"(ATS job prefix stripping produced duplicates).",
                           file=sys.stderr )
                    return 1
                results[ key ] = data
    parent = os.path.dirname( os.path.abspath( args.output ) )
    if parent:
        os.makedirs( parent, exist_ok=True )
    with open( args.output, "w" ) as f:
        json.dump( results, f, indent=2, sort_keys=True )
    scope = " `_iterative` geosx" if args.iterative else ""
    print( f"Harvested iteration counts from {len( results )}{scope} log file(s) under {root}" )
    if not results:
        print( "Error: no per-solve iteration lines found; "
               "are the tests running with LinearSolverParameters logLevel >= 1?",
               file=sys.stderr )
        return 1
    return 0


def firstSequenceDiff( first, second ):
    """Return (index, a, b) for the first mismatch, or None if sequences match."""
    n = min( len( first ), len( second ) )
    for i in range( n ):
        if first[ i ] != second[ i ]:
            return i, first[ i ], second[ i ]
    if len( first ) != len( second ):
        a = first[ n ] if n < len( first ) else None
        b = second[ n ] if n < len( second ) else None
        return n, a, b
    return None


def compareExact( first, second, key ):
    marks = []
    a, b = first[ key ], second[ key ]
    seqA, seqB = a[ "iterations" ], b[ "iterations" ]
    if a[ "numFailedSolves" ] != b[ "numFailedSolves" ]:
        marks.append( f"failed-solve counts differ "
                      f"({a[ 'numFailedSolves' ]} vs {b[ 'numFailedSolves' ]})" )
    diff = firstSequenceDiff( seqA, seqB )
    if diff is not None:
        index, valA, valB = diff
        extra = abs( len( seqA ) - len( seqB ) )
        extraNote = f" (length {len( seqA )} vs {len( seqB )})" if extra else ""
        marks.append( f"solve {index + 1}: {valA} vs {valB}{extraNote}" )
        remaining = 0
        start = index + 1
        n = min( len( seqA ), len( seqB ) )
        for i in range( start, n ):
            if seqA[ i ] != seqB[ i ]:
                remaining += 1
        remaining += extra
        if remaining:
            marks.append( f"and {remaining} later difference(s)" )
    return marks


def compareSlack( first, second, key, args ):
    marks = []
    a, b = first[ key ], second[ key ]
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
    return marks


def compare( args ):
    with open( args.first ) as f:
        first = json.load( f )
    with open( args.second ) as f:
        second = json.load( f )

    common = sorted( set( first ) & set( second ) )
    onlyFirst = sorted( set( first ) - set( second ) )
    onlySecond = sorted( set( second ) - set( first ) )

    if not first or not second:
        print( "Error: one or both harvests are empty." )
        return 1

    if not common:
        print( "Error: no common log files between the two harvests." )
        return 1

    violations = []
    header = f"{'log file':60s} {'solves':>13s} {'total its':>15s} {'max its':>13s}"
    print( header )
    print( "-" * len( header ) )
    for key in common:
        a, b = first[ key ], second[ key ]
        if args.exact_sequence:
            marks = compareExact( first, second, key )
        else:
            marks = compareSlack( first, second, key, args )

        status = "FAIL: " + "; ".join( marks ) if marks else "ok"
        print( f"{key:60s} {a[ 'numSolves' ]:>6d}/{b[ 'numSolves' ]:<6d} "
               f"{a[ 'totalIterations' ]:>7d}/{b[ 'totalIterations' ]:<7d} "
               f"{a[ 'maxIterations' ]:>6d}/{b[ 'maxIterations' ]:<6d}  {status}" )
        if marks:
            violations.append( key )

    if onlyFirst or onlySecond:
        print( f"\nNote: {len( onlyFirst )} log file(s) only in {args.first}, "
               f"{len( onlySecond )} only in {args.second}." )
        if args.require_identical_keys:
            for key in onlyFirst:
                print( f"  only in {args.first}: {key}" )
            for key in onlySecond:
                print( f"  only in {args.second}: {key}" )
            print( "Error: harvest key sets differ." )
            return 1

    print( f"\nCompared {len( common )} log file(s): "
           f"{len( common ) - len( violations )} ok, {len( violations )} violation(s)." )
    return 1 if violations else 0


class CompareLinearSolverIterationsTests( unittest.TestCase ):
    def testCanonicalKeyStripsAtsPrefix( self ):
        key = canonicalKey(
            "poromechanics/PoroElastic_Terzaghi_iterative_smoke_01/"
            "0135.geosx_PoroElastic_Terzaghi_iterative_smoke_01_1_geos_.log",
            True )
        self.assertEqual(
            key,
            "poromechanics/PoroElastic_Terzaghi_iterative_smoke_01/"
            "geosx_PoroElastic_Terzaghi_iterative_smoke_01_1_geos_.log" )

    def testShouldHarvestIterativeGeosxOnly( self ):
        class Args:
            iterative = True
            geosx_only = True
        rel = "poromechanics/PoroElastic_Terzaghi_iterative_smoke_01/" \
              "0135.geosx_PoroElastic_Terzaghi_iterative_smoke_01_1_geos_.log"
        self.assertTrue( shouldHarvest( rel, os.path.basename( rel ), Args ) )
        restart = rel.replace( "geosx_", "python3_" ).replace( "_geos_.log",
                                                               "_restartcheck_.log" )
        self.assertFalse( shouldHarvest( restart, os.path.basename( restart ), Args ) )
        other = "poromechanics/PoroElastic_Terzaghi_direct_01/0001.geosx_foo.log"
        self.assertFalse( shouldHarvest( other, os.path.basename( other ), Args ) )

    def testFirstSequenceDiff( self ):
        self.assertIsNone( firstSequenceDiff( [ 1, 2 ], [ 1, 2 ] ) )
        self.assertEqual( firstSequenceDiff( [ 1, 2, 3 ], [ 1, 9, 3 ] ), ( 1, 2, 9 ) )
        self.assertEqual( firstSequenceDiff( [ 1, 2 ], [ 1, 2, 3 ] ), ( 2, None, 3 ) )

    def testExactCompareDetectsPerSolveDrift( self ):
        first = { "t": { "iterations": [ 8, 8, 9 ],
                         "numSolves": 3,
                         "totalIterations": 25,
                         "maxIterations": 9,
                         "numFailedSolves": 0 } }
        second = { "t": { "iterations": [ 8, 9, 9 ],
                          "numSolves": 3,
                          "totalIterations": 26,
                          "maxIterations": 9,
                          "numFailedSolves": 0 } }
        marks = compareExact( first, second, "t" )
        self.assertTrue( any( "solve 2: 8 vs 9" in m for m in marks ) )


def selftest( _args ):
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(
        CompareLinearSolverIterationsTests )
    result = unittest.TextTestRunner( verbosity=2 ).run( suite )
    if result.wasSuccessful():
        with tempfile.TemporaryDirectory() as tmp:
            logDir = os.path.join( tmp, "poromechanics", "Poro_iterative_01" )
            os.makedirs( logDir )
            logPath = os.path.join( logDir, "0135.geosx_Poro_iterative_01_1_geos_.log" )
            with open( logPath, "w" ) as f:
                f.write( "Linear Solver | Success | Unknowns: 10 | Iterations: 12 | Final rel. res.: 1e-9\n" )
            out = os.path.join( tmp, "out.json" )
            class HarvestArgs:
                directory = tmp
                output = out
                iterative = True
                geosx_only = True
                strip_ats_prefix = True
            if harvest( HarvestArgs ) != 0:
                return 1
            with open( out ) as f:
                data = json.load( f )
            expected = "poromechanics/Poro_iterative_01/geosx_Poro_iterative_01_1_geos_.log"
            if expected not in data or data[ expected ][ "iterations" ] != [ 12 ]:
                print( f"selftest harvest mismatch: {data}", file=sys.stderr )
                return 1
        return 0
    return 1


def addHarvestFlags( parser ):
    parser.add_argument( "--iterative", action="store_true",
                         help="restrict to `_iterative` geosx logs, strip ATS job "
                              "number prefixes, and fail if nothing is harvested" )
    parser.add_argument( "--geosx-only", action="store_true",
                         help="ignore logs whose filename does not contain 'geosx'" )
    parser.add_argument( "--strip-ats-prefix", action="store_true",
                         help="strip leading `<job>.` from harvested log names" )


def main():
    parser = argparse.ArgumentParser( description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter )
    sub = parser.add_subparsers( dest="command", required=True )

    h = sub.add_parser( "harvest", help="collect iteration counts from log files" )
    h.add_argument( "directory", help="directory to scan recursively" )
    h.add_argument( "-o", "--output", required=True, help="output JSON file" )
    addHarvestFlags( h )
    h.set_defaults( func=harvest )

    c = sub.add_parser( "compare", help="compare two harvested JSON files" )
    c.add_argument( "first" )
    c.add_argument( "second" )
    c.add_argument( "--exact-sequence", action="store_true",
                    help="require identical per-solve iteration sequences (zero slack)" )
    c.add_argument( "--require-identical-keys", action="store_true",
                    help="fail if the two harvests do not contain the same log keys" )
    c.add_argument( "--total-rel-tol", type=float, default=0.1,
                    help="(slack mode) allowed relative difference in total iterations "
                         "per log (default 0.1)" )
    c.add_argument( "--total-abs-tol", type=int, default=5,
                    help="(slack mode) absolute total-iteration slack below which the "
                         "relative criterion is not applied (default 5)" )
    c.add_argument( "--max-solve-abs-tol", type=int, default=5,
                    help="(slack mode) allowed difference in the worst per-solve "
                         "iteration count (default 5)" )
    c.set_defaults( func=compare )

    t = sub.add_parser( "selftest", help="run built-in unit tests" )
    t.set_defaults( func=selftest )

    args = parser.parse_args()
    if getattr( args, "iterative", False ):
        args.geosx_only = True
        args.strip_ats_prefix = True
    if getattr( args, "exact_sequence", False ):
        args.require_identical_keys = True
    return args.func( args )


if __name__ == "__main__":
    sys.exit( main() )

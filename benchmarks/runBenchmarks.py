import os
import glob
import xml.etree.ElementTree as ElementTree
import time
import datetime
import functools
import shutil
import subprocess
import sys
import argparse


class Status:
    """ Represents the status of a job. """

    NOT_SUBMITTED = "NOT_SUBMITTED"
    SUBMITTED = "SUBMITTED"
    SUCCESS = "SUCCESS"
    FAILURE = "FAILURE"


def getMostCubeLikeRepresentation( n ):
    """
    Of all possible boxes of integer dimensions and volume n return the dimensions of the most cube-like box.

    Args:
        n: The volume of the box.

    Returns:
        The dimensions of the most cube-like box returned as a triple.
    """
    minTriple = None
    minSurfaceArea = 1e99
    for i in range( 1, n + 1 ):
        if n % i != 0:
            continue

        n_over_i = n / i
        for j in range( 1, n_over_i + 1 ):
            if n_over_i % j != 0:
                continue
            
            k = n_over_i / j
            surfaceArea = 2 * i * j + 2 * i * k + 2 * j * k
            if surfaceArea < minSurfaceArea:
                minSurfaceArea = surfaceArea
                minTriple = ( i, j, k )

    return minTriple


def parseListFromString( listString ):
    """
    Given a comma separated string construct a list. The list may optionally be enclosed in {}. 

    Args:
        listString: The string to parse.

    Returns:
        A list of strings.
    """
    newList = listString.strip(" {}").split(",")
    return map( lambda x: x.strip(), newList )


def createDirectory( dirPath, clean=False ):
    """
    Create a new directory structure, it is ok if the directory already exists.

    Args:
        dirPath: Path to the directory to create.
        clean: True if the directory is to be cleaned.
    """
    if os.path.exists( dirPath ):
        if clean:
            shutil.rmtree( dirPath )
        elif os.path.isdir( dirPath ):
            return
        else:
            raise Exception( "{} already exists but it is not a directory.".format( dirPath ) )

    os.makedirs( dirPath )


class Machine:
    """
    An class that implements an abstract computing system.

    Concrete classes are intended to derive from Machine and provide
    specific submission commands for that machine.

    Attributes:
        name: The name of the Machine.
    """

    def __init__( self, name ):
        """
        Initialize a Machine with the given name.

        Args:
            self: The Machine to initialize.
            name: The name of the Machine.
        """
        self.name = name

    def getSumbissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        This method needs to be overridden in the derived class.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes.
        
        Returns:
            A list containing the submission commands.
        """
        raise NotImplementedError( "This method must be overridden in the derived class." )

    def printProgress( self ):
        """
        Print the status of the submission queue.

        This method needs to be overridden in the derived class.

        Args:
            self: The Machine to get the status of.
        """
        raise NotImplementedError( "This method must be overridden in the derived class." )

    def hasCudaGPU( self ):
        """ Return True if this Machine has a CUDA GPU available."""
        return False


class SlurmMachine( Machine ):
    """ A class that implements a machine that uses Slurm. """

    def getSumbissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes, may be None if not specified.
        
        Returns:
            A list containing the submission commands.
        """
        command = ["srun", "-N", nodes, "-n", tasks]
        if threadsPerTask is not None:
            command += ["-c", threadsPerTask]
        if timeLimit is not None:
            command += ["-t", timeLimit]

        return command

    def printProgress( self ):
        """
        Print the status of the submission queue.

        Args:
            self: The Machine to get the status of.
        """
        print( subprocess.check_output( [ "squeue", "-u", os.environ[ "USER" ] ], stderr=subprocess.STDOUT ) )


class LallocMachine( Machine ):
    """ A class that implements a machine that uses Lalloc. """

    def getSumbissionCommand( self, nodes, tasks, threadsPerTask, timeLimit ):
        """
        Return a list containing the command necessary to launch a job with the given configuration on this machine.

        Args:
            self: The Machine to get the submission command for.
            nodes: The number of nodes the job will use.
            tasks: The number of tasks in the job.
            threadsPerTask: The threads per task, may be None if not specified.
            timeLimit: The time limit of the job in minutes, may be None if not specified.
        
        Returns:
            A list containing the submission commands.
        """
        # TODO(corbett5): something about the bsub "-Is" argument messes up the printing while jobs are running.
        command = ["lalloc", nodes, "lrun", "-n", tasks, "-g", 1]
        if threadsPerTask is not None:
            command += ["--threadsPerTask", threadsPerTask ]
        if timeLimit is not None:
            command += ["-W", timeLimit]
        
        return command

    def printProgress( self ):
        """
        Print the status of the submission queue.

        Args:
            self: The Machine to get the status of.
        """
        print( subprocess.check_output( [ "bjobs" ], stderr=subprocess.STDOUT ) )
    
    def hasCudaGPU( self ):
        """ Return True."""
        return True


def getMachine():
    """ Return a Machine object corresponding to the current machine. """
    hostName = os.environ[ "HOSTNAME" ]
    if hostName.startswith( "quartz" ):
        return SlurmMachine( "quartz" )
    elif hostName.startswith( "lassen" ):
        return LallocMachine( "lassen" )
    else:
        raise EnvironmentError( "Cannot determine the current machine." )


class Benchmark:
    """
    An class that implements a single benchmark run of GEOSX.

    Attributes:
        geosxPath: The path to the GEOSX executable to run.
        xmlPath: The path to the input XML file to pass to GEOSX.
        name: The name of the benchmark.
        nodes: The number of nodes to run the benchmark with.
        tasks: The number of tasks to run the benchmark with.
        threadsPerTask: The number of threads per task to run the benchmark with.
            May be None if not specified.
        timeLimit: The time limit to give the scheduler in minutes, or None if not specified.
        outputDir: The directory where the benchmark is run.
        outputFile: The path to the file containing the standard output and
            standard error from the benchmark.
        runCommand: A list of arguments which appended to the submission command provides
            the full command for running this benchmark.
        process: The subproccess associated with the benchmark.
        status: The status of the benchmark.
    """

    def __init__( self, outputDir, geosxPath, xmlPath, name, nodes, tasks, threadsPerTask, timeLimit, args, autoPartition ):
        """
        Initialize a Benchmark.

        Arguments:
            self: The Benchmark to initialize.
            outputDir: The top level directory where the benchmarks are run.
            geosxPath: The path to the GEOSX executable to run.
            xmlPath: The path to the input XML file to pass to GEOSX.
            name: The name of the benchmark.
            nodes: The number of nodes to run the benchmark with.
            tasks: The number of tasks to run the benchmark with.
            threadsPerTask: The number of threads per task to run the benchmark with.
                May be None if not specified.
            timeLimit: The time limit to give the scheduler in minutes, or None if if not specified.
            args: List of extra arguments to pass to GEOSX.
            autoPartition: If true then partition arguments are generated and passed to GEOSX.
        """
        self.geosxPath = os.path.abspath( geosxPath )
        self.xmlPath = os.path.abspath( xmlPath )
        self.name = name
        self.nodes = nodes
        self.tasks = tasks
        self.threadsPerTask = threadsPerTask
        self.timeLimit = timeLimit

        xmlName = os.path.splitext( os.path.basename( self.xmlPath ) )[ 0 ]
        self.outputDir = os.path.abspath( os.path.join( outputDir, xmlName, "{}_{}".format( self.name, self.nodes ) ) )

        self.outputFile = os.path.join( self.outputDir, "output.txt" )

        self.runCommand = [self.geosxPath, "-n", "{}/{}".format( xmlName, self.name ), "-i", self.xmlPath]

        self.runCommand += args

        if autoPartition:
            pX, pY, pZ = getMostCubeLikeRepresentation( tasks )
            self.runCommand += ["-x", pX, "-y", pY, "-z", pZ]

        self.process = None
        self.status = Status.NOT_SUBMITTED

    def submit( self, machine ):
        """
        Submit to be run on the given machine.

        Arguments:
            self: The Benchmark to submit.
            machine: The machine to run on.
        """
        createDirectory( self.outputDir, True )

        submissionCommand = machine.getSumbissionCommand( self.nodes, self.tasks, self.threadsPerTask, self.timeLimit )

        submissionCommand += self.runCommand

        if machine.hasCudaGPU():
            submissionCommand += ["-t", "spot,profile.mpi,profile.cuda"]
        else:
            submissionCommand += ["-t", "spot,profile.mpi"]
        
        submissionCommand = map( str, submissionCommand )
        with open( self.outputFile, "w" ) as outputFile:
            outputFile.write( "{}\n\n".format( " ".join( submissionCommand ) ) )
            self.process = subprocess.Popen( submissionCommand, cwd=self.outputDir, stdout=outputFile, stderr=subprocess.STDOUT )

        self.status = Status.SUBMITTED

        print( "Submitted {}".format( self ) )
    
    def getOutputFile( self ):
        """
        Return the path to the file containing the standard output and error of the Benchmark.

        Arguments:
            self: The Benchmark to get the output of.
        """
        return self.outputFile

    def getTimingFile( self ):
        """
        Return the path to the timing file produced by the Benchmark.

        Arguments:
            self: The Benchmark to get the timing file of.
        """
        caliFiles = glob.glob( os.path.join( self.outputDir, "*.cali" ) )
        if caliFiles:
            if len( caliFiles ) > 1:
                raise Exception( "The directory {} contains multiple .cali files!".format( self.outputDir ) )
            return caliFiles[ 0 ]
        else:
            return None

    def hasCompleted( self ):
        """
        Return True iff the Benchmark has completed.

        Arguments:
            self: The Benchmark to check.
        """
        jobExited = self.status == Status.SUBMITTED and self.process.poll() is not None
        if self.getTimingFile() is not None:
            self.status = Status.SUCCESS
            print( "Completed {}".format( self ) )
        elif jobExited:
            self.status = Status.FAILURE
            print( "Failed {}".format( self ) )

        return self.status == Status.SUCCESS or self.status == Status.FAILURE

    def kill( self ):
        """
        Kill the Benchmark.

        Arguments:
            self: The Benchmark to kill.
        """
        self.process.kill()
        self.process.wait()
        self.hasCompleted()

    def __str__( self ):
        """
        Return a string represtation of the Benchmark.

        Arguments:
            self: The Benchmark to stringify.
        """
        return "Benchmark( OutputDir: {}, nodes: {}, tasks: {}, threadsPerTask: {}, timeLimit {} )".format( self.outputDir, self.nodes, self.tasks, self.threadsPerTask, self.timeLimit )


def printProgress( machine, benchmarks, startTime ):
    """
    Print the status of the submitted benchmarks.

    Arguments:
        machine: The Machine the benchmarks were submitted on.
        benchmarks: A list of the submitted Benchmarks.
        startTime: The time of the benchmarks were submitted at.
    """
    print("")

    runTime = int( time.time() - startTime )
    print( "Running for {}".format( datetime.timedelta( seconds=runTime ) ) )

    statusCount = { Status.NOT_SUBMITTED: 0, Status.SUBMITTED: 0, Status.SUCCESS: 0, Status.FAILURE: 0 }
    for benchmark in benchmarks:
        statusCount[ benchmark.status ] += 1

    for status, count in statusCount.iteritems():
        print( "{}: {}".format( status, count ) )

    machine.printProgress()
    sys.stdout.flush()
    print("")


def allCompleted( benchmarks ):
    """
    Return True if all the benchmarks have completed.

    Arguments:
        benchmarks: A list of the submitted Benchmarks.
    """
    return functools.reduce( lambda result, benchmark : result and benchmark.hasCompleted(), benchmarks, True )


def submitAllAndWait( machine, benchmarks, timeLimit ):
    """
    Submit all of the benchmarks and wait until the complete or the time limit expires.

    Arguments:
        machine: The Machine to submit the benchmarks to.
        benchmarks: A list of the Benchmarks to submit.
        timeLimit: The number of minutes to wait after submitting the benchmarks.
    """
    startTime = time.time()
    for benchmark in benchmarks:
        benchmark.submit( machine )

    while not allCompleted( benchmarks ):
        for _ in range( 60 ):
            time.sleep( 1 )
            if allCompleted( benchmarks ):
                break

        printProgress( machine, benchmarks, startTime )

        minutesRun = ( time.time() - startTime ) / 60
        if minutesRun > timeLimit:
            print( "The time limit of {} minutes has been exceeded, killing all unfinished jobs.".format( timeLimit ) )
            for benchmark in benchmarks:
                benchmark.kill()


def getBenchmarksFromXML( xmlFilePath, machine, outputDir, geosxPath ):
    """
    Return a list of benchmarks created for the current Machine from the given XML file.

    Arguments:
        xmlFilePath: The path to the XML file.
        machine: The Machine to run the benchmarks on.
        outputDir: The top level directory all the benchmarks are to be run in.
        geosxPath: The path to the GEOSX executable to run.
    """
    benchmarks = []

    tree = ElementTree.parse( xmlFilePath )
    matchingElements = tree.findall( "./Benchmarks/{}/Run".format( machine.name ) )

    for elem in matchingElements:
        name = elem.get( "name" )
        if name is None:
            raise Exception( "Expected the XML element to have an attribute 'name'." )

        nodes = elem.get( "nodes" )
        if nodes is None:
            raise Exception( "Expected the XML element to have an attribute 'nodes'." )
        nodes = int( nodes )

        tasksPerNode = elem.get( "tasksPerNode" )
        if tasksPerNode is None:
            raise Exception( "Expected the XML element to have an attribute 'tasksPerNode'." )
        tasksPerNode = int( tasksPerNode )

        threadsPerTask = elem.get( "threadsPerTask" )
        if threadsPerTask is not None:
            threadsPerTask = int( threadsPerTask )

        timeLimit = elem.get( "timeLimit" )
        if timeLimit is not None:
            timeLimit = int( timeLimit )

        args = elem.get( "args" )
        if args is None:
            args = []
        else:
            args = args.split()

        autoPartition = elem.get( "autoPartition" )
        if autoPartition is None:
            autoPartition = False
        elif autoPartition in ("Off", "On"):
            autoPartition = autoPartition == "On"
        else:
            raise Exception( "Option for autoPartition not recognized: {}. Valid options are 'On' and 'Off'.".format( autoPartition ) )

        strongScaling = elem.get( "strongScaling" )
        if strongScaling is None:
            benchmarks.append( Benchmark( outputDir, geosxPath, xmlFilePath, name, nodes, nodes * tasksPerNode, threadsPerTask, timeLimit, args, autoPartition ) )
        else:
            strongScaling = parseListFromString( strongScaling )
            for scale in strongScaling:
                totalNodes = nodes * int( scale )
                totalTasks = totalNodes * tasksPerNode
                benchmarks.append( Benchmark( outputDir, geosxPath, xmlFilePath, name, totalNodes, totalTasks, threadsPerTask, timeLimit, args, autoPartition ) )

    return benchmarks


def getBenchmarksFromDirectory( benchmarkDir, machine, outputDir, geosxPath ):
    """
    Return a list of benchmarks created for the current Machine from the XML files in the given directory.

    Arguments:
        benchmarkDir: The directory containing the XML files to run.
        machine: The Machine to run the benchmarks on.
        outputDir: The top level directory all the benchmarks are to be run in.
        geosxPath: The path to the GEOSX executable to run.
    """
    benchmarks = []
    for fileName in os.listdir( benchmarkDir ):
        if fileName.endswith( ".xml" ):
            benchmarks += getBenchmarksFromXML( os.path.join( benchmarkDir, fileName ), machine, outputDir, geosxPath )

    return benchmarks


def main():
    """ Parse the command line arguments, submit all found benchmarks and wait for them to finish. """

    benchmarkDir = os.path.dirname( os.path.realpath( __file__ ) )

    timeLimit = 60
    parser = argparse.ArgumentParser()
    parser.add_argument( "geosxPath", help="The path to the GEOSX executable to benchmark.")
    parser.add_argument( "outputDirectory", help="The parent directory to run the benchmarks in.")
    parser.add_argument( "-t", "--timeLimit", type=int, help="Time limit for the entire script in minutes, the default is {}.". format( timeLimit ), default=timeLimit )
    parser.add_argument( "-o", "--timingCollectionDir", help="Directory to copy the timing files to." )
    parser.add_argument( "-e", "--errorCollectionDir", help="Directory to copy the output from any failed runs to." )
    args = parser.parse_args()

    geosxPath = os.path.abspath( args.geosxPath )
    if not os.path.isfile( geosxPath ) or not os.access( geosxPath, os.X_OK):
        raise ValueError( "geosxPath is not an executable!" )

    outputDir = os.path.abspath( args.outputDirectory )
    
    timeLimit = args.timeLimit
    
    timingCollectionDir = args.timingCollectionDir
    if timingCollectionDir is not None:
        timingCollectionDir = os.path.abspath( timingCollectionDir )

    errorCollectionDir = args.errorCollectionDir
    if errorCollectionDir is not None:
        errorCollectionDir = os.path.abspath( errorCollectionDir )

    machine = getMachine()

    benchmarks = getBenchmarksFromDirectory( benchmarkDir, machine, outputDir, geosxPath )

    print( "Benchmarking GEOSX found at {}".format( geosxPath ) )
    print( "Results will be written to {}".format( outputDir ) )
    if timingCollectionDir is not None:
        print( "Timing files will be added to {}".format( timingCollectionDir ) )
    if errorCollectionDir is not None:
        print( "Output from failed benchmarks will be put in {}".format( errorCollectionDir ) )
    print( "The time limit is {} minutes.".format( timeLimit) )
    print( "Found {} benchmarks.".format( len( benchmarks ) ) )
    print( "" )

    submitAllAndWait( machine, benchmarks, timeLimit )

    # Copy the timing files from successful benchmarks to a new directory if asked.
    if timingCollectionDir is not None:
        createDirectory( timingCollectionDir )
        for benchmark in benchmarks:
            if benchmark.status == Status.SUCCESS:
                shutil.copy2( benchmark.getTimingFile(), timingCollectionDir )

    # Copy the output from the failed benchmarks to a new directory if asked.
    if errorCollectionDir is not None:
        createDirectory( errorCollectionDir, True )
        for benchmark in benchmarks:
            if benchmark.status != Status.SUCCESS and os.path.isfile( benchmark.getOutputFile() ):
                outputFile = benchmark.getOutputFile()
                newLocation = os.path.join( errorCollectionDir, os.path.relpath( outputFile, outputDir ) )
                shutil.copy2( outputFile, newLocation )
    
    success = True
    for benchmark in benchmarks:
        success = success and benchmark.status == Status.SUCCESS
    
    return not success


if __name__ == "__main__" and not sys.flags.interactive:
    sys.exit(main())

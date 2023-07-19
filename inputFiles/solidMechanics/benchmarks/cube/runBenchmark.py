import subprocess
import re


def cgThroughput( executable, inputFile, numRuns ):

  print( "cgThroughput( ", executable, ", ", inputFile, ", ", numRuns, " )" )
  throughput = []
  for i in range(numRuns):
    print( "run ", i, " of ", numRuns )
    proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', executable, '-i', inputFile], stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'Throughput'], stdin=proc1.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, err = proc2.communicate()
    throughput.append( float( out.decode('utf-8').split(' ')[1] ) )

  #print( throughput )
  return max(throughput)


def femKernelThroughput_nvprof( executable, inputFile, numRuns ):
  print( "femKernelThroughput_nvprof( ", executable, ", ", inputFile, ", ", numRuns, " )" )

  throughput = []
  for i in range(numRuns):
    print( "run ", i, " of ", numRuns )
    proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', 'nvprof', executable, '-i', inputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'Small'], stdin=proc1.stderr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, err = proc2.communicate()

    result = re.split('([mu])', out.decode('utf-8').split()[7] )
    unit = result[1]
    if unit=='u':
      value = float(result[0]) * 1.0e-6;
    elif unit=='m':
      value = float(result[0]) * 1.0e-3;
    elif unit=='s':
      value = float(result[0]) * 1.0e-3;    
    else:
      raise Exception("invalid unit")

    throughput.append( value )

  return min(throughput)




numDofs = [ 11*11*11*3,
            101*11*11*3,
            101*101*11*3,
            101*101*101*3,
            1001*101*101*3,
            1001*1001*101*3,
            1001*1001*1001*3 
            ]

cgThroughputs = []
cgThroughputs.append( cgThroughput( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_111.xml', 10 ) )
cgThroughputs.append( cgThroughput( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_211.xml', 10 ) )
cgThroughputs.append( cgThroughput( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_221.xml', 10 ) )
cgThroughputs.append( cgThroughput( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_222.xml', 10 ) )
cgThroughputs.append( cgThroughput( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_322.xml', 10 ) )


kernelThroughputs = []
kernelThroughputs.append( femKernelThroughput_nvprof( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_111.xml', 10 ) )
kernelThroughputs.append( femKernelThroughput_nvprof( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_211.xml', 10 ) )
kernelThroughputs.append( femKernelThroughput_nvprof( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_221.xml', 10 ) )
kernelThroughputs.append( femKernelThroughput_nvprof( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_222.xml', 10 ) )
kernelThroughputs.append( femKernelThroughput_nvprof( 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_322.xml', 10 ) )



for case in range(len(cgThroughputs)):
  print( "{0:10d} {1:6.2f} {2:6.2f}".format( numDofs[case], cgThroughputs[case], numDofs[case]/kernelThroughputs[case]/1.0e6 ) )
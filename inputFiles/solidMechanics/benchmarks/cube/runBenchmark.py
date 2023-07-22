import subprocess
import re


def cgThroughput( executable, inputFile, numRuns ):

  #print( "cgThroughput( ", executable, ", ", inputFile, ", ", numRuns, " )" )
  throughput = []
  for i in range(numRuns):
    #print( "run ", i, " of ", numRuns )
    #proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', executable, '-i', inputFile], stdout=subprocess.PIPE)
    proc1 = subprocess.Popen([executable, '-i', inputFile], stdout=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'MaxThroughput'], stdin=proc1.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, err = proc2.communicate()
    throughput.append( float( out.decode('utf-8').split(' ')[1] ) )

  #print( throughput )
  return max(throughput)


def femKernelTime_nvprof( executable, inputFile, numRuns ):
  #print( "femKernelTime_nvprof( ", executable, ", ", inputFile, ", ", numRuns, " )" )

  kernelTime = []
  for i in range(numRuns):
  #  print( "run ", i, " of ", numRuns )
    #proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', 'nvprof', executable, '-i', inputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1 = subprocess.Popen(['nvprof', executable, '-i', inputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
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
      value = float(result[0]);
    else:
      raise Exception("invalid unit")

    kernelTime.append( value )

  return min(kernelTime)


def femKernelTime_caliper( executable, inputFile, numRuns ):

  kernelTime = []
  for i in range(numRuns):
    #proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', 'nvprof', executable, '-i', inputFile, '-t runtime-report,max_column_width=400' ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1 = subprocess.Popen([executable, '-i', inputFile, '-t', 'runtime-report'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'kernelLaunch'], stdin=proc1.stderr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #out, err = proc1.communicate()
    out2, err2 = proc2.communicate()
    value = float(out2.decode('utf-8').split()[1]) / 1000.0 ;

    kernelTime.append( value )

  return min(kernelTime)






numDofs = [ 
            11*11*11*3,
            101*11*11*3,
            101*101*11*3,
            101*101*101*3,
            1001*101*101*3,
            1001*501*101*3,
            1001*1001*101*3,
            1001*1001*251*3 
            ]


runList = [
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_111.xml', 5],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_211.xml', 5],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_221.xml', 5],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_222.xml', 5],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_322.xml', 5],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_322plus.xml', 5 ],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_332.xml', 5 ],
            [ 'bin/geosx', '../inputFiles/solidMechanics/benchmarks/cube/cube_332plus.xml', 5 ]
          ]


cgThroughputs = []
print( "CG Throughput")
print( "   #Dofs     MDof/s")
for case in range( 0, len(runList) ):
  cgThroughputs.append( cgThroughput( runList[case][0], runList[case][1], runList[case][2] ) )
  print( "{0:10d} {1:6.2f}".format( numDofs[case], cgThroughputs[case] ) )


kernelTimes = []
print( "Kernel Throughput")
print( "   #Dofs     MDof/s")
for case in range( 0, len(runList) ):
  kernelTimes.append( femKernelTime_caliper( runList[case][0], runList[case][1], runList[case][2] ) )
  print( "{0:10d} {1:6.2f}".format( numDofs[case], numDofs[case]/kernelTimes[case]/1.0e6 ) )


# for case in range(len(kernelTimes)):
#   print( "{0:10d} {1:6.2f} {2:6.2f}".format( numDofs[case], numDofs[case]/kernelTimes[case]/1.0e6 ) )
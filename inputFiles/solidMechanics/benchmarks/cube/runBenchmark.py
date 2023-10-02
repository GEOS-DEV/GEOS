import subprocess
import re
import os
import os.path
import time

# on perlmutter this may be required to disable default counters...then reenable them.
#  subprocess.Popen("srun --ntasks-per-node 1 dcgmi profile --pause", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#  subprocess.Popen("srun --ntasks-per-node 1 dcgmi profile --resume", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def getLaunchCommands():
  nerscSystemName = os.environ.get('NERSC_HOST')
  lcSystemName = os.environ.get('LCSCHEDCLUSTER')
  systemName = os.environ.get('HOSTNAME')
  noEnvResult = os.environ.get('thisEnvironmentVariableDoesNotExist')

  if nerscSystemName != noEnvResult and lcSystemName != noEnvResult:
    print( 'both nerscSystemName (',nerscSystemName,') and lcSystemName (',lcSystemName,') are set. This does not make sense' )
  else:
    if nerscSystemName != noEnvResult:
      systemName = nerscSystemName

    if lcSystemName != noEnvResult:
      systemName = lcSystemName


  launchCommand = []

  if systemName == 'lassen':
    launchCommand = [ 'lrun', '-N1', '-n1', '-g1' ]
  # elif systemName == 'perlmutter':
  #   launchCommand = [ 'srun', '-N1', '-n1', '-G1', ]

  return launchCommand



def cgThroughput( executable, inputFile, numRuns ):
  launchCommand = getLaunchCommands()
  throughput = []
  for _ in range(numRuns):
    #print( "run ", i, " of ", numRuns )
    proc1 = subprocess.Popen( launchCommand + [executable, '-i', inputFile], stdout=subprocess.PIPE )
    proc2 = subprocess.Popen( ['grep', 'MaxThroughput'], stdin=proc1.stdout,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
#    time.sleep(5.0)
    out, _ = proc2.communicate()
    throughput.append( float( out.decode('utf-8').split(' ')[1] ) )

  #print( throughput )
  return max(throughput)


def femKernelTime_nvprof( executable, inputFile, numRuns ):
  #print( "femKernelTime_nvprof( ", executable, ", ", inputFile, ", ", numRuns, " )" )

  kernelTime = []
  for _ in range(numRuns):
  #  print( "run ", i, " of ", numRuns )
    #proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', 'nvprof', executable, '-i', inputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1 = subprocess.Popen(['nvprof', executable, '-i', inputFile], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'Small'], stdin=proc1.stderr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, _ = proc2.communicate()

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


def femKernelTime_nsys( executable, inputFile, numRuns ):
  #print( "femKernelTime_nsys( ", executable, ", ", inputFile, ", ", numRuns, " )" )

  version_cmd = ['nsys', '--version']
  proc0 = subprocess.Popen(version_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  out, _ = proc0.communicate()
  versionOutput = out.decode('utf-8');
  #print( versionOutput )
  if '2022' in versionOutput:
    reportScript = 'gpukernsum'
  elif '2023' in versionOutput:
    reportScript = 'cuda_gpu_kern_sum'
  else:
    reportScript = 'cuda_gpu_kern_sum'
    #print('nsys version is not 2022 or 2023. Assuming same 2023 syntax!!')

  launchCommand = getLaunchCommands()

  kernelTime = []
  for _ in range(numRuns):
    report_fname = "/tmp/" + os.path.basename(inputFile) + ".nsys-rep" # TODO make robust
    profile_cmd = launchCommand + ['nsys', 'profile', '-o', report_fname, '-f', 'true', executable, '-i', inputFile]
    stats_cmd = ['nsys', 'stats', '--force-overwrite', 'true', '--force-export', 'true', '--timeunit', 'nsec', '--report', reportScript, '--format', 'column:nohdr:nolimit', report_fname]

    subprocess.run(profile_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).check_returncode()
    proc1 = subprocess.Popen(stats_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    time.sleep(1.0)
    proc2 = subprocess.Popen(['grep', 'SmallStrainResidual'], stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, _ = proc2.communicate()
    stringValue = out.decode('utf-8').split()[5];
    value = float( stringValue.replace(',','') ) * 1.0e-9
    kernelTime.append( value )

  return min(kernelTime)

def femRoofline_ncu( executable, inputFile, numRuns ):
  launchCommand = getLaunchCommands()

  kernelTime = []
  for _ in range(numRuns):
    report_fname = "/tmp/" + os.path.basename(inputFile) + ".nsys-rep" # TODO make robust
    profile_cmd = launchCommand + ['ncu',
                                   '--kernel-name-base', 'demangled',
                                   '--kernel-id', '::regex:SmallStrainResidual:',
                                   '--target-processes', 'all',
                                   '--metrics', 'dram__bytes_write.sum,dram__bytes_read.sum,smsp__sass_thread_inst_executed_op_dadd_pred_on.sum,smsp__sass_thread_inst_executed_op_dmul_pred_on.sum,smsp__sass_thread_inst_executed_op_dfma_pred_on.sum,gpu__time_duration.min',
                                   '--launch-count', '10',
                                   '--print-summary','per-gpu',
                                   '--print-units','base',
                                    executable, '-i', inputFile]


    proc1 = subprocess.Popen(profile_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', '-A', '7', 'Metric Name'], stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1.stdout.close() # Allow proc1 to receive a SIGPIPE if proc2 exits.
    out, _ = proc2.communicate()
    stringValue = out.decode('utf-8')
#    print(stringValue)
    lines = stringValue.split('\n')

    memory = 0
    dflops = 0
    duration = 0

    for c in range(2,len(lines)):
      if lines[c]:
        line = lines[c].split()
#        print(line)
        name = line[0]
        units = line[1]
        minVal = float(line[2])
        if any([x in name for x in ['dram__bytes_read.sum','dram__bytes_write.sum']]):
          memory += minVal
        elif any([x in name for x in ['smsp__sass_thread_inst_executed_op_dadd_pred_on.sum','smsp__sass_thread_inst_executed_op_dmul_pred_on.sum']]):
          dflops += minVal
        elif 'smsp__sass_thread_inst_executed_op_dfma_pred_on.sum' in name:
          dflops += 2*minVal
        elif 'duration' in name:
          duration = minVal

    return memory, dflops, duration
    # value = float( stringValue.replace(',','') ) * 1.0e-9
    # kernelTime.append( value )

  return min(kernelTime)

def femRoofline_rocprof( executable, inputFile, numRuns ):
  import pandas as pd
  launchCommand = getLaunchCommands()
  rocprof_infile =  os.path.join(os.path.dirname(os.path.realpath(__file__)),'rocprof-input.txt')
  rocprof_outfile = os.path.join(os.path.dirname(os.path.realpath(__file__)),'rocprof-output.csv')

  profile_cmd = launchCommand + ['rocprof',
                                  '-i', rocprof_infile,
                                  '-o', rocprof_outfile,
                                  '--timestamp', 'on',
                                  executable, '-i', inputFile]
  subprocess.check_output(profile_cmd)
  df = pd.read_csv(rocprof_outfile)
  hmb_memory = 32 * ( df['TCC_EA_RDREQ_32B_sum'] + df['TCC_EA_WRREQ_sum'] - df['TCC_EA_WRREQ_64B_sum'] ) + \
               64 * ( df['TCC_EA_RDREQ_sum'] - df['TCC_EA_RDREQ_32B_sum'] + df['TCC_EA_WRREQ_64B_sum'] )

  # total_f32_flop = 64 * ( df['SQ_INSTS_VALU_ADD_F32'] + df['SQ_INSTS_VALU_MUL_F32'] + df['SQ_INSTS_VALU_TRANS_F32'] +  2 * df['SQ_INSTS_VALU_FMA_F32'] ) + 512 * ( df['SQ_INSTS_VALU_MFMA_MOPS_F32'] )
  total_f64_flop = 64 * ( df['SQ_INSTS_VALU_ADD_F64'] + df['SQ_INSTS_VALU_MUL_F64'] + df['SQ_INSTS_VALU_TRANS_F64'] +  2 * df['SQ_INSTS_VALU_FMA_F64'] ) + 512 * ( df['SQ_INSTS_VALU_MFMA_MOPS_F64'] )

  return hmb_memory.mean(), total_f64_flop.mean(), femKernelTime_rocprof( executable, inputFile, numRuns )

def femKernelTime_rocprof( executable, inputFile, numRuns ):
  kernelTime = []
  for _ in range(numRuns):
    profile_cmd = ['rocprof', '--hip-trace', executable, '-i', inputFile ]

    subprocess.run(profile_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).check_returncode()

    import json
    with open('results.json') as fin:
      results = json.loads( fin.read() )
    rtimes = []
    for result in results['traceEvents']:
      if 'name' in result.keys():
        if 'SmallStrainResidual' in result['name']:
          rtimes.append( int( result['args']['DurationNs'] ) )

    kernelTime.append( min( rtimes ) )

  return min(kernelTime)

def femKernelTime_caliper( executable, inputFile, numRuns ):

  kernelTime = []
  for _ in range(numRuns):
    #proc1 = subprocess.Popen(['lrun', '-N1', '-n1', '-g1', 'nvprof', executable, '-i', inputFile, '-t runtime-report,max_column_width=400' ], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc1 = subprocess.Popen([executable, '-i', inputFile, '-t', 'runtime-report'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc2 = subprocess.Popen(['grep', 'kernelLaunch'], stdin=proc1.stderr, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    #out, err = proc1.communicate()
    out2, _ = proc2.communicate()
    value = float(out2.decode('utf-8').split()[1]) / 1000.0

    kernelTime.append( value )

  return min(kernelTime)




benchmark_dir = os.path.dirname(os.path.realpath(__file__))

runList = [
            #  ( 'bin/geosx', 'cube_111.xml',     2, 11*11*11*3 ), // can't use this for rocm, we "converge" in ~840 steps +=3 and rocprof can't merge metrics from runs with different number of kernel launches
             ( 'bin/geosx', 'cube_211.xml',     2, 101*11*11*3 ),
             ( 'bin/geosx', 'cube_221.xml',     2, 101*101*11*3 ),
             ( 'bin/geosx', 'cube_222.xml',     2, 101*101*101*3 ),
             ( 'bin/geosx', 'cube_322.xml',     2, 1001*101*101*3 ),
             ( 'bin/geosx', 'cube_322plus.xml', 2, 1001*501*101*3 )
            # ( 'bin/geosx', 'cube_332.xml',     2, 1001*1001*101*3 ),
            # ( 'bin/geosx', 'cube_332plus.xml', 2, 1001*1001*251*3 )
          ]

import sys

if 'cg' in sys.argv:
  print( "CG Throughput")
  print( "   #Dofs     MDof/s")
  for executable, inputFile, numRuns, numDofs in runList:
    throughput = cgThroughput( executable, os.path.join(benchmark_dir, inputFile), numRuns )
    print( "{0:10d} {1:>8.2f}".format( numDofs, throughput ) )

if 'caliper' in sys.argv:
  print( "Kernel Throughput (caliper)")
  print( "   #Dofs     MDof/s")
  for executable, inputFile, numRuns, numDofs in runList:
    kernelTime = femKernelTime_caliper( executable, os.path.join(benchmark_dir, inputFile), numRuns )
    print( "{0:10d} {1:>8.2f}".format( numDofs, numDofs/kernelTime/1.0e6 ) )

if 'rocprof' in sys.argv:
  print( "Kernel Throughput (rocprof)")
  print( "   #Dofs     MDof/s")
  for executable, inputFile, numRuns, numDofs in runList:
    kernelTime = femKernelTime_rocprof( executable, os.path.join(benchmark_dir, inputFile), numRuns ) / 1e9
    print( "{0:10d} {1:>8.2f}".format( numDofs, numDofs/kernelTime/1.0e6 ) )

if 'nsys' in sys.argv:
  print( "Kernel Throughput (nsys)")
  print( "   #Dofs     MDof/s")
  for executable, inputFile, numRuns, numDofs in runList:
    kernelTime = femKernelTime_nsys( executable, os.path.join(benchmark_dir, inputFile), numRuns )
    print( "{0:10d} {1:>8.2f}".format( numDofs, numDofs/kernelTime/1.0e6 ) )

if 'ncu' in sys.argv:
  print( "Roofline")
  print( "     #Dofs HBM(bytes)       FLOP     dur(s)         AI    TFLOP/s")

  for executable, inputFile, numRuns, numDofs in runList:
    memory, dflop, duration = femRoofline_ncu( executable, os.path.join(benchmark_dir, inputFile), 1 )
    print( "{0:10d} {1:>8.4e} {2:>8.4e} {3:>8.4e} {4:>8.4e} {5:>8.4e}".format( numDofs, memory, dflop, duration/1e9, dflop/memory, dflop/(duration/1e9)/1e12 ) )

if 'rocprof-roofline' in sys.argv:
  print( "Roofline")
  print( "     #Dofs HBM(bytes)       FLOP     dur(s)         AI    TFLOP/s")

  for executable, inputFile, numRuns, numDofs in runList:
    memory, dflop, duration = femRoofline_rocprof( executable, os.path.join(benchmark_dir, inputFile), 1 )
    print( "{0:10d} {1:>8.4e} {2:>8.4e} {3:>8.4e} {4:>8.4e} {5:>8.4e}".format( numDofs, memory, dflop, duration/1e9, dflop/memory, dflop/(duration/1e9)/1e12 ) )
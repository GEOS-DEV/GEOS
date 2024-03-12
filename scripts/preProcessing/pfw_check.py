# -*- coding: utf-8 -*-
from __future__ import print_function    # (at top of module)
from __future__ import division
from __future__ import unicode_literals
import datetime                               # used for date stamp
import os                                     # operating sys commands, pwd, etc.
import subprocess                             # lets you call msub and get jobid etc
import sys                                    # to access command arguments.
import importlib                              # to import userDefs

def reverse_readline(filename, buf_size=8192):
    """A generator that returns the lines of a file in reverse order"""
    with open(filename, 'rb') as fh:
        segment = None
        offset = 0
        fh.seek(0, os.SEEK_END)
        file_size = remaining_size = fh.tell()
        while remaining_size > 0:
            offset = min(file_size, offset + buf_size)
            fh.seek(file_size - offset)
            buffer = fh.read(min(remaining_size, buf_size)).decode(encoding='utf-8')
            remaining_size -= buf_size
            lines = buffer.split('\n')
            # The first line of the buffer is probably not a complete line so
            # we'll save it and append it to the last line of the next buffer
            # we read
            if segment is not None:
                # If the previous chunk starts right from the beginning of line
                # do not concat the segment to the last line of new chunk.
                # Instead, yield the segment first 
                if buffer[-1] != '\n':
                    lines[-1] += segment
                else:
                    yield segment
            segment = lines[0]
            for index in range(len(lines) - 1, 0, -1):
                if lines[index]:
                    yield lines[index]
        # Don't yield None if the file was empty
        if segment is not None:
            yield segment

# ##################################################################
# # READ IN JOB FILE TO GET USER PARAMETERS FOR THE SIMULATION
# ##################################################################
username = os.getenv('LOGNAME') or os.getenv('USER') or os.getenv('LNAME') or os.getenv('USERNAME')
userDefsFile = str('userDefs_'+str(username))

try:
  f = open(userDefsFile+".py")
  # Do something with the file
except IOError:
  print("Please create a version of userDefs_OUN.py consistent with your user name and geos location")
finally:
  f.close()

userDefs = importlib.import_module(userDefsFile)
geosPath = userDefs.geosPath

# read timeStamp as argument.  If it's missing give a warning.
if len(sys.argv) != 3:
  print("usage : pfw_check.py inputFile jobID \nYou must specify the timeStamp as argument")
  sys.exit(1)

# slurm file will be slurm-jobID.out
inputFile = sys.argv[1]
jobID = sys.argv[2]

#import inputFile as job    # Input parameters in separate file.
job = importlib.import_module(inputFile)

pfw = job.pfw

# parameters for batch job:
mBatch = pfw["mBatch"]
mBank = pfw["mBank"]
mWallTime = pfw["mWallTime"]
mNodes = pfw["mNodes"]
mCores = pfw["mCores"]
xpar = pfw["xpar"]
ypar = pfw["ypar"]
zpar = pfw["zpar"]
stopTime = pfw["stopTime"]

[wH,wM,wS]=mWallTime.split(":")
mWallTimeMinutes=int(wH)*60+int(wM)
mWallTime = mWallTimeMinutes
partition = "pbatch"
if pfw["runDebug"]:
  mWallTimeMinutes=min(mWallTimeMinutes,60) # in minutes
  partition="pdebug"

# Current working directory
PWD = os.getcwd()

needsRestart = False
jobExitedForUnknownReason = True

# lastTimestep = 0.0
slurmFile = "slurm-"+str(jobID)+".out"
for line in reverse_readline(slurmFile):
  if 'Job complete' in line:
    jobExitedForUnknownReason = False
    break

  if 'Job exited early' in line or 'TIME LIMIT' in line:
    jobExitedForUnknownReason = False
    needsRestart= True
    break

if jobExitedForUnknownReason:
  print("Job interrupted for unknown reason.")
  sys.exit()

geosInputFileName = 'mpm_'+inputFile.replace('pfw_input_',"")+'.xml'

# file name prefix is date stamped
now = datetime.datetime.now()
timeStamp = str(now.year)[-2:]+str(now.month).zfill(2)+str(now.day).zfill(2) + \
            str(now.hour).zfill(2)+str(now.minute).zfill(2)+str(now.second).zfill(2)

if needsRestart:
  # check for an 8 digit restart file.
  nq = 16
  restartFile = ""
  while (restartFile == "" and nq > 5):
    restartFile = subprocess.Popen('ls -dtr *_restart_'+'?'*nq+' | tail -1',shell=True,stdout=subprocess.PIPE).communicate()[0]
    if ( len(restartFile) == 0 ):
      nq = nq - 1
      restartFile = ""
    elif ( str(restartFile, 'UTF8')[-5]=="."):
      print("skipping ",str(restartFile, 'UTF8'))
      nq = nq - 1
      restartFile = ""

  restartFile = str(restartFile, 'UTF8').replace('.root','')
  print('last restartFile = ',restartFile)

  if restartFile=="restart_000000":
    print('Simulation ended before restart file written, reduce restart interval')
    sys.exit(1)


  # ===========================================
  #  Make a SLURM script to run GEOS.
  # ===========================================

  # This will run the job.
  slurmScript = """#!/bin/bash
#SBATCH -t """+str(mWallTimeMinutes)+"""
#SBATCH -N """+str(mNodes)+"""
#SBATCH -p """+ partition +"""
#SBATCH -A """+mBank+"""

srun -n """+str(mCores)+""" """+geosPath+""" -i """+geosInputFileName+""" -x """ + str(xpar) + """ -y """ + str(ypar) + """ -z """ + str(zpar) + """ -r """+restartFile+"""
"""

  fileName = timeStamp+"_restartGEOS.sh"
  file = open(fileName, 'w')
  file.write(slurmScript)
  file.close()

  output = subprocess.Popen(["sbatch", fileName], stdout=subprocess.PIPE).communicate()[0]
  output = str(output, 'UTF8')
  output = output.strip('Submitted batch job ')
  jobID = output.strip()
  print('run_check output = ',output.strip())

  slurmScript = """#!/bin/bash
#SBATCH -t 00:02:00
#SBATCH -N 1
#SBATCH -p """+ partition +"""
#SBATCH -A """+mBank+"""
#SBATCH --dependency=afterany:"""+jobID+"""

python3 pfw_check.py """+inputFile+""" """+jobID+"""
"""

  fileName = timeStamp+"_runCheck.sh"
  file = open(fileName, 'w')
  file.write(slurmScript)
  file.close()
  
  subprocess.call(["sbatch",fileName], cwd=PWD)


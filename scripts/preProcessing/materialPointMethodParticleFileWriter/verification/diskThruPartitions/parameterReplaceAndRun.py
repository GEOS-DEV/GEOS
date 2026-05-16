# -*- coding: utf-8 -*-
from __future__ import print_function    # (at top of module)
from __future__ import division
from __future__ import unicode_literals
import os                                     # operating sys commands, pwd, etc.
from subprocess import call                   # lets you call shell commands, i.e. call(["ln", "-s",".","run_dir"])
import subprocess                             # lets you call msub and get jobid

# ==============================================================
# 1D compaction with continuum porosity, relaxation and recovery
# ==============================================================
pfwDirectory="/g/g19/homel1/pfwx/"
templateFileName="pfw_input_diskThruSpace_rXXXX.py"
templateFileLocation='/usr/workspace/homel1/pfw_geosx_temp/verification/diskThruPartitions'
replaceStrings=["XXXX"]

import platform
lassen = 'lassen' in platform.node()
if lassen:
  machineName = 'lassen' 
  runCommand = 'bsub'
else:
  machineName = 'quartz'
  runCommand = 'sbatch'

values = [2,4,8,16,32,64]

for v0 in values:
	newFileName = templateFileName.replace(replaceStrings[0],str(v0).replace(".","p"))
	print(newFileName)
	with open(templateFileName,"rt") as fin:
		with open(newFileName,"wt") as fout:
			for line in fin:
				new_line = line.replace(replaceStrings[0],str(v0))
				fout.write(new_line)
	with open(pfwDirectory+"runClean_"+machineName+"_Template.sh","rt") as fin:
		with open(pfwDirectory+"runClean_"+machineName+"_Temp.sh","wt") as fout:
			for line in fin:
				new_line = line.replace("REPLACEME",(newFileName.replace("pfw_input_","")).replace(".py",""))
				new_line = new_line.replace("FILELOCATION",templateFileLocation)
				fout.write(new_line)
	# copy the input file to a time-stamped file
	call([runCommand,"runClean_"+machineName+"_Temp.sh"],cwd=pfwDirectory)


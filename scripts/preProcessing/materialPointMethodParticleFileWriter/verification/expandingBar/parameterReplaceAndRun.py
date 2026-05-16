# -*- coding: utf-8 -*-
from __future__ import print_function    # (at top of module)
from __future__ import division
from __future__ import unicode_literals
import os                                     # operating sys commands, pwd, etc.
from subprocess import call                   # lets you call shell commands, i.e. call(["ln", "-s",".","run_dir"])
import subprocess                             # lets you call msub and get jobid

# ==============================================================
# Expanding Ring study on strength, fracture energy, and strain rate
# ==============================================================
#templateFileName="pfw_input_expandingRing_chiumenti_sigmaXXXX_GfYYYY_v0ZZZZ.py"
templateFileName="pfw_input_expandingBar_chiumenti_Weibull_sigmaXXXX_GfYYYY_v0ZZZZ.py"
#templateFileName="pfw_input_expandingRing_chiumenti_sigmaXXXX_GfYYYY_v0ZZZZ_r56.py"
templateFileLocation="/g/g19/homel1/workspace/particleFileWriter_exportControlled/particleFileWriter/verification/expandingBar"
replaceStrings=["XXXX","YYYY","ZZZZ"]

# initialRadialVelocity = 0.01592
# density = 2.75						# mg/mm^3
# bulk = 275./( 3.*( 1. - 2.*0.3) )	# GPa
# shear = 275./( 2.*( 1. + 0.3) )		# GPa
# failureStrength = 0.30				# GPa
# energyReleaseRate = 1.0e-5 			# mg/us^2
#values = [[0.3000],[0.00001],[0.0016000,0.0028,0.0050,0.0089,0.0160,0.0284,0.0505,0.0899,0.1600,0.2845,0.5059,0.8997,1.6000]]

values = [[0.3000],[0.00001],[0.1,1.0]]
#values = [[0.3000],[0.00001],[0.0160]]

for v0 in values[0]:
	for v1 in values[1]:
		for v2 in values[2]:
			newFileName = ((templateFileName.replace(replaceStrings[0],str(v0).replace(".","p"))).replace(replaceStrings[1],str(v1).replace(".","p"))).replace(replaceStrings[2],str(v2).replace(".","p"))
			print(newFileName)

			with open(templateFileName,"rt") as fin:
				with open(newFileName,"wt") as fout:
					for line in fin:
						new_line = line.replace(replaceStrings[0],str(v0))
						new_line = new_line.replace(replaceStrings[1],str(v1))
						new_line = new_line.replace(replaceStrings[2],str(v2))
						fout.write(new_line)

			with open("/g/g19/homel1/workspace/particleFileWriter_exportControlled/particleFileWriter/runClean_ruby_Template.sh","rt") as fin:
				with open("/g/g19/homel1/workspace/particleFileWriter_exportControlled/particleFileWriter/runClean_homel1_Temp.sh","wt") as fout:
					for line in fin:
						new_line = line.replace("REPLACEME",(newFileName.replace("pfw_input_","")).replace(".py",""))
						new_line = new_line.replace("FILELOCATION",templateFileLocation)
						fout.write(new_line)
			# copy the input file to a time-stamped file
			call(["sbatch","runClean_homel1_Temp.sh"],cwd="/g/g19/homel1/workspace/particleFileWriter_exportControlled/particleFileWriter/")

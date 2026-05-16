import numpy as np   
import sys
sys.path.insert(1,'/g/g19/homel1/pfwx/particleFileWriter')
import pfw_geometryObjects as geom              
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os
import re
import importlib
from cycler import cycler
import logging
from pfw_analysis import *

# INPUTS ------------------------------------------------------------------------------------------
# This script will search for jobs with names specified in paths variable in the locations specified here
runLocations=['/p/lustre1/homel1/geosxRuns/']

# Can specify any number of jobs by their name (pfw_input_<NAME>.py)
paths = [ 'geoModel_HYDcreep'
        ]

plotReactions = True # If false will only output results to console window
readBoxSums = True #True # Read the boxAverage csv file output by jbos

# filter data?
filterByMedian = True
windowSize = 5

# MISC SCRIPT VARIABLES -------------------------------------------------------------------------------

plt.rcParams.update({'font.size': 22})
legend_font_size = 8

# BEGIN EXECUTION -------------------------------------------------------------------------------------

job_dirs = format_file_paths(runLocations, paths)

print('Found following job files:')
for d in job_dirs:
  print('\t' + d)

# Keep track of which files were read correctly and which could not be read due to error
read_files = []
errored_files = []

# Curves from different jobs will have unique color
color_index = np.linspace(0, 1, len(job_dirs))
colors = [ cm.rainbow(x) for x in color_index ]

fig, axes = plt.subplots(3, 1,figsize=(12,8))

for i,file in enumerate(job_dirs):
  print()
  print('Processing', file)
  print()

  jobObj = MPMJob(file)
  filename = jobObj.job_name

  try:
    jobObj.read_reaction_file()
    jobObj.read_from_box_average_file()
    jobObj.applyPostProcess("all", RemoveNonMonotonicEntries(jobObj.fields["Time"].getData()))
    if (filterByMedian):
      jobObj.applyPostProcess("all", MedianFilter(windowSize))
       
    plt.figure(1)
    if readBoxSums:
      # Plot domain stress
      Btime = jobObj.fields["BTime"].getData()
      Bsxx = jobObj.fields["BSxx"].getData()
      Bsyy = jobObj.fields["BSyy"].getData()
      Bszz = jobObj.fields["BSzz"].getData()
      Bp = (-1/3)*(Bsxx+Bsyy+Bszz)    

      axes[0].plot(Btime,Bsxx,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bsxx")
      axes[0].plot(Btime,Bsyy,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bsyy")
      axes[0].plot(Btime,Bszz,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bszz")
      axes[0].plot(Btime,Bp,linestyle='-',color=colors[i],linewidth=3,alpha=0.5,label=filename+"_Bp")

      # Plot domain plastic strain
      Bepxx = jobObj.fields["Bepxx"].getData()
      Bepyy = jobObj.fields["Bepyy"].getData()
      Bepzz = jobObj.fields["Bepzz"].getData()
      Bevp = (Bepxx+Bepyy+Bepzz)

      axes[2].plot(Btime,Bepxx,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepxx")
      axes[2].plot(Btime,Bepyy,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepyy")
      axes[2].plot(Btime,Bepzz,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepzz")
      axes[2].plot(Btime,Bevp,linestyle='-',color=colors[i],linewidth=3,alpha=0.5,label=filename+"_Bevp")

    if plotReactions:
      engineeringStress = False
      jobObj.compute_domain_stress(engineeringStress)
      time = jobObj.fields["Time"].getData()
      
      # pot domain deformation
      F00 = jobObj.fields["F00"].getData()
      F11 = jobObj.fields["F11"].getData()
      F22 = jobObj.fields["F22"].getData()     
      axes[1].plot(time,F00,linestyle='-',color=colors[i],linewidth=1,label=filename+"_F00")
      axes[1].plot(time,F11,linestyle='--',color=colors[i],linewidth=1,label=filename+"_F11")
      axes[1].plot(time,F22,linestyle=':',color=colors[i],linewidth=1,label=filename+"_F22")
      
      # plot domain stress computed from reactions:
      Rsxx = jobObj.fields["Rsxx"].getData()
      Rsyy = jobObj.fields["Rsyy"].getData()
      Rszz = jobObj.fields["Rszz"].getData()
      axes[0].plot(time,Rsxx,linestyle='-',color=colors[i],linewidth=1,label=filename+"_Rsxx")
      axes[0].plot(time,Rsyy,linestyle='-',color=colors[i],linewidth=1,label=filename+"_Rsyy")
      axes[0].plot(time,Rszz,linestyle='-',color=colors[i],linewidth=1,label=filename+"_Rszz")
      
    read_files.append(file)

  except BaseException:
    logging.exception("Could not read file: " + file)
    print()
    errored_files.append(file)

# Format plots:
plt.figure(1)
axes[0].set_ylabel(r'Stress')
axes[0].tick_params(axis='both', which='major')
axes[0].tick_params(axis='both', which='minor')
axes[0].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[0].grid()

axes[1].set_ylabel(r'F')
axes[1].tick_params(axis='both', which='major')
axes[1].tick_params(axis='both', which='minor')
axes[1].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[1].grid()

axes[2].set_xlabel('Time (μs)')
axes[2].set_ylabel(r'ep')
axes[2].tick_params(axis='both', which='major')
axes[2].tick_params(axis='both', which='minor')
axes[2].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[2].grid()

fig.tight_layout()
fig.savefig("geoHYDcreep.png", bbox_inches="tight")

print()
print("******************** Done processing jobs ********************")
print()
if len(read_files) != 0:
  print("Read jobs:")
  for file in read_files:
    print('\t' + file)
print()
if len(errored_files) != 0: 
  print("Could not read files:")
  for file in errored_files:
    print('\t' + file)
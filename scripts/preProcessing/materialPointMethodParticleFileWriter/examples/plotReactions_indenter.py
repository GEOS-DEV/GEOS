import numpy as np   
import pfw_geometryObjects as geom              
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
import os
import re
import importlib
from cycler import cycler
import logging

from pfw_analysis import *

# INPUTS ------------------------------------------------------------------------------------------
# This script will search for jobs with names specified in paths variable in the locations specified here
runLocations=['/p/lustre1/crook5/geosxRuns/', 
              '/p/lustre2/crook5/geosxRuns/',
              '/p/lustre2/crook5/geosxDev/']

# Can specify any number of jobs by their name (pfw_input_<NAME>.py)
paths = [   
        # 'spherical_indenter_2D',
        # 'spherical_indenter_2D_implicit',
        # 'spherical_indenter_2D_implicit_ppcy5',
        # 'spherical_indenter_2D_implicit_g144',
        # 'spherical_indenter_2D_elasticOnly',

        # 'sphr2DIndenter_elastic_explicit_g120',
        # 'sphr2DIndenter_elastic_implicit_g120',


        # 'sphr2DIndenter_elastic_explicit_g60',

	      # 'sphrIndenter2D_g60_frictionless_explicit',
	      # 'sphrIndenter2D_g60_frictionless_implicit',
	      # 'sphrIndenter2D_g60_slip_explicit',
	      # 'sphrIndenter2D_g60_slip_implicit',
	      # 'sphrIndenter2D_g60_stick_explicit',
	      # 'sphrIndenter2D_g60_stick_implicit',

        # 'sphrIndenter2D_g100_frictionless_explicit',
	      # 'sphrIndenter2D_g100_frictionless_implicit',
	      # 'sphrIndenter2D_g100_slip_explicit',
	      # 'sphrIndenter2D_g100_slip_implicit',
	      # 'sphrIndenter2D_g100_stick_explicit',
	      # 'sphrIndenter2D_g100_stick_implicit',

        'sphrIndenter2D_g120_frictionless_explicit',
	      'sphrIndenter2D_g120_frictionless_implicit',
	      # 'sphrIndenter2D_g120_slip_explicit',
	      # 'sphrIndenter2D_g120_slip_implicit',
	      # 'sphrIndenter2D_g120_stick_explicit',
	      # 'sphrIndenter2D_g120_stick_implicit',
        ]

plotReactions = True # If false will only output results to console window
readBoxSums = False #True # Read the boxAverage csv file output by jbos
livePlot = False # Leave plotting window open, otherwise just saves image

writePostprocessFile = True # writes data to file

useEngineeringStrain=True

flipCompression = False #True # Flips stress for compression simulations
filterStresses = False # Removes spike artifacts in reaction data above a stress threshold
stressMaxThreshold = np.inf # GPa
stressMinThreshold = -stressMaxThreshold

enforceAxesLimits = True # Handy for manually specifying axes limits in plots
yStressMax=700
yStressMin=-5

filterByMedian = True
windowSize = 5

subSampleData = True
numSubSamples = 2000

# END INPUTS --------------------------------------------------------------------------------------


# MISC SCRIPT VARIABLES -------------------------------------------------------------------------------

plt.rcParams.update({'font.size': 22})

labels = []
volumes = []
modulus = []
youngsModulus = []
poissonRatio = []
bulkModulus = []

maxS1 = []
maxVM = []
strengths = []
loadings = []
strains = []
failed = []
particles = []

legend_font_size = 8
maxStress = 0
minStress = 0

deformationGradMin = 1.0
deformationGradMax = 1.0

# END MISC SCRIPT VARIABLES ---------------------------------------------------------------------------


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

fig, axes = plt.subplots(2, 1,figsize=(12,8))
axes = np.append(axes, axes[0].twinx())

for i,file in enumerate(job_dirs):
  print()
  print('Processing', file)
  print()

  jobObj = MPMJob(file)
  filename = jobObj.job_name

  try:
    jobObj.read_reaction_file()
    jobObj.applyPostProcess("all", RemoveNonMonotonicEntries(jobObj.fields["Time"].getData()))
    jobObj.applyPostProcess("all", MedianFilter(windowSize))
    # jobObj.applyPostProcess("all", SubSample(method="nearest", numSamples=numSubSamples))
       
    if plotReactions:
      plt.figure(1)

      time = jobObj.fields["Time"].getData()
      displacement = -1000.0*jobObj.domainHeight0*(jobObj.fields["F11"].getData()-1.0)
      load = -1000.0*jobObj.fields["Ryp"].getData()/jobObj.domainLength0

      # Plot domain deformation gradient
      axes[0].plot(time,displacement,linestyle='-',color=colors[i],linewidth=1,label=filename+"_u")
      axes[2].plot(time,load,linestyle='--',color=colors[i],linewidth=1, label=filename+"_P")

      axes[1].plot(displacement, load,linestyle='-',color=colors[i],linewidth=1, label=filename)

    read_files.append(file)

  except BaseException:
    logging.exception("Could not read file: " + file)
    print()
    errored_files.append(file)

# # Plot Analytical form  
# K_d = 54.0
# G_d = 44.0
# K_s = 62.0
# G_s = 25.0

# R = 0.001
# b = 0.1

# E_d = 9*K_d*G_d/(3*K_d+G_d)
# E_s = 9*K_s*G_s/(3*K_s+G_s)
# v_d = (3*K_d-2*G_d)/(2*(3*K_d+G_d))
# v_s = (3*K_s-2*G_s)/(2*(3*K_s+G_s))
# E_reduced = 1/((1-v_s**2)/E_s+(1-v_d**2)/E_d)

# F = np.linspace(0.001,0.0025,100)
# L = 4*F*R/(np.pi*b**2*E_reduced)

# u = 1e3*(F / (np.pi*E_reduced*L) )*( np.log( 1e-6*4*np.pi*E_reduced*R*L/F )-1)
# print(F, u)
# axes[1].plot(u, F)


plt.figure(1)
axes[0].set_xlabel('Time (μs)')
axes[0].set_ylabel(r'u (μm)')
axes[0].tick_params(axis='both', which='major')
axes[0].tick_params(axis='both', which='minor')
axes[2].set_ylabel(r'P (mN)')
axes[0].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))

axes[1].set_xlabel(r'u (μm)')
axes[1].set_ylabel(r'P (mN)')
# axes[1].set_xlim([0, 0.01])
# axes[1].set_ylim([0, 0.25])
axes[1].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[1].grid()

fig.tight_layout()
fig.savefig("indenter_reactions.png", bbox_inches="tight")

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

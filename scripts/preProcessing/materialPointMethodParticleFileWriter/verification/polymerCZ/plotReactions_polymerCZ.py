import numpy as np            
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
sys.path.insert(0,'/usr/workspace/crook5/pfw_dev/')
import os
import re
import importlib
from cycler import cycler
import logging

from pfw_analysis import *
import pfw_geometryObjects as geom    



# INPUTS ------------------------------------------------------------------------------------------
# This script will search for jobs with names specified in paths variable in the locations specified here
runLocations=['/p/lustre1/crook5/geosxDev/', 
              '/p/lustre2/crook5/geosxDev/']

# Can specify any number of jobs by their name (pfw_input_<NAME>.py)
paths = [              
        'polymerCZ',
        ]

plotReactions = True # If false will only output results to console window
useEngineeringStrain=True

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

    # Remove nonmonotonic entries (happens sometimes when job restarts)
    monotonicFilter = RemoveNonMonotonicEntries(jobObj.fields["Time"].getData())
    jobObj.applyPostProcess("all", monotonicFilter)
    if filterByMedian:
        print("Performing median filter...")
        medianFilter = MedianFilter(windowSize)
        jobObj.applyPostProcess("all", medianFilter)
    
    if subSampleData:
        print("Performing sub sample...")
        subSample = SubSample(method="nearest", numSamples=numSubSamples)
        jobObj.applyPostProcess("all", subSample)

    jobObj.compute_domain_stress()
    jobObj.compute_domain_strain(engineeringStrain=True)

    if plotReactions:
        plt.figure(1)

        # experimental data for 1e-3/s and 0C, strain is log strain, stress is MPa
        dataTension_strain=np.array([0.,0.0152041,0.0257091,0.0399285,0.0613519,0.0899947,0.125734,0.190025,0.254301,0.404224,0.586254,0.768264,0.960942,1.06796,1.06986])
        dataTension_stress=np.array([0.1,17.1782,25.7759,29.3165,29.6106,25.7825,24.5412,24.3057,24.7979,28.7834,35.8106,46.3524,65.5892,83.2086,2.09566])

        dataCompression_strain=np.array([0,0.00456397,0.012768,0.0249448,0.0471481,0.0683025,0.0874232,0.126676,0.165935,0.21224,0.326003,0.386414,0.44078,0.500182,0.500155,0.479909,0.462677,0.450462,0.446298])
        dataCompression_stress=np.array([0.1,3.98692,11.183,21.6455,31.7168,34.7841,33.196,31.3032,30.5627,29.8357,30.8236,32.6245,33.3545,34.8972,28.9904,13.4961,5.99795,2.37367,0.9284])

        axes[1].plot(-1.0*dataCompression_strain*jobObj.sampleHeight,-1.0*dataCompression_stress,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Data - Compression')
        axes[1].plot(dataTension_strain*jobObj.sampleHeight,dataTension_stress,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Data - Tension')



        time = jobObj.fields["Time"].getData()
        eyy = jobObj.fields["eyy"].getData()
        uyy = jobObj.sampleHeight*eyy
        F00 = jobObj.fields["F00"].getData()
        F11 = jobObj.fields["F11"].getData()
        F22 = jobObj.fields["F22"].getData()
        rsxm = 1000.0*jobObj.fields["Rsxm"].getData()
        rsxp = 1000.0*jobObj.fields["Rsxp"].getData()
        rsym = 1000.0*jobObj.fields["Rsym"].getData()
        rsyp = 1000.0*jobObj.fields["Rsyp"].getData()
        rsxx = 1000.0*jobObj.fields["Rsxx"].getData()
        rsyy = 1000.0*jobObj.fields["Rsyy"].getData()

        # Plot domain deformation gradient
        axes[0].plot(time, F00,linestyle='-',color=colors[i],linewidth=1,label=filename+"_F00")
        axes[0].plot(time, F11,linestyle='--',color=colors[i],linewidth=1, label=filename+"_F11")
        axes[0].plot(time, F22,linestyle='-.',color=colors[i],linewidth=1, label=filename+"_F22")

        axes[2].plot(time, rsxm,linestyle='-',color=colors[i],linewidth=1,label=filename+"_rsxx_-x")
        axes[2].plot(time, rsxp,linestyle='--',color=colors[i],linewidth=1, label=filename+"_rsxx_+x")
        axes[2].plot(time, rsym,linestyle='-',color=lighten_color(colors[i], 0.5),linewidth=1,label=filename+"_rsyy_-y")
        axes[2].plot(time, rsyp,linestyle='--',color=lighten_color(colors[i], 0.5),linewidth=1, label=filename+"_rsyy_+y")
        axes[2].plot(time, rsxx,linestyle='-',color=lighten_color(colors[i], 1.0),linewidth=1,label=filename+"_rsxx")
        axes[2].plot(time, rsyy,linestyle='-',color=lighten_color(colors[i], 0.5),linewidth=1,label=filename+"_rsyy")

        axes[1].plot(uyy, rsyy,linestyle='-',color=colors[i],linewidth=1, label=filename+"_rsyy")
        axes[1].plot(uyy, rsym,linestyle='--',color=colors[i],linewidth=1, label=filename+"_rsym")
        axes[1].plot(uyy, rsyp,linestyle='-.',color=colors[i],linewidth=1, label=filename+"_rsyp")

    read_files.append(file)

  except BaseException:
    logging.exception("Could not read file: " + file)
    print()
    errored_files.append(file)

plt.figure(1)
axes[0].set_xlabel('Time (us)')
axes[0].set_ylabel(r'F (-)')
axes[0].tick_params(axis='both', which='major')
axes[0].tick_params(axis='both', which='minor')

axes[2].set_ylabel(r'$\sigma$ (MPa)')

axes[1].set_xlabel(r'Displacement (mm)')
axes[1].set_ylabel(r'$\sigma$ (MPa)')
axes[1].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[1].grid()

fig.tight_layout()
fig.savefig("polymerCZ_reactions.png", bbox_inches="tight")

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


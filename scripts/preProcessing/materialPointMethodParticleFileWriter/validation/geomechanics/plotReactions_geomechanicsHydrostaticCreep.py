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
#for interpolating the time and strain data
from scipy import interpolate

# INPUTS ------------------------------------------------------------------------------------------
# This script will search for jobs with names specified in paths variable in the locations specified here
runLocations=['/p/lustre1/homel1/geosxRuns/']

# Can specify any number of jobs by their name (pfw_input_<NAME>.py)
paths = [ 'geomechanicsHydrostaticCreep_HP03'
        ]

shorthand='HP03_IS03'

A = 0.222    # creep: phi_e = A*exp(-p/B), note that A is 0-pressure equilbrium porosity and phi_i = 1 -exp(-p3) = 0.225
B = 0.00084     # creep: phi_e = A*exp(-p/B), B is GPa
equilibriumPorosityPressureExponent = 2.15

p3 = 0.255

phi_i = 1-np.exp(-p3)
print('phi_i=',phi_i)

stopTime = 1000.  #simulation end time
physicalStopTime= 90.78*24*60*60*1000000 # for micro seconds

timeScaleToHrs = 90.78/1000  # this will be used to scale creep-rate term.


#initialize plots
##plot strain versus time data from the experiment
timeData=[]
strainData=[]

with open ('strain_vs_time_'+shorthand+'.csv') as csvfile:
	for line in csvfile:
		line = line.strip();
		values = line.split(',');
		
		#for some reason, there are other characters in the first column. Strip them.
		#values[0]=re.sub(r'[^0-9.]','',values[0]);

		values[0]=float(values[0]);
		values[1]=float(values[1]);
		timeData.append(values[0]);
		strainData.append(values[1]);



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

#fig, axes = plt.subplots(4, 1,figsize=(24,30))
fig, axes= plt.subplots(4, 1,figsize=(16,16))

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
      
      BMatVol = jobObj.fields["BMatVol"].getData()
      BJ = BMatVol / BMatVol[0]
      Bev = np.log(BJ)

      # axes[0].plot(Btime,-Bsxx,linestyle=':',color=colors[i],linewidth=1,label='-'+filename+"_Bsxx")
      # axes[0].plot(Btime,-Bsyy,linestyle=':',color=colors[i],linewidth=1,label='-'+filename+"_Bsyy")
      # axes[0].plot(Btime,-Bszz,linestyle=':',color=colors[i],linewidth=1,label='-'+filename+"_Bszz")
      axes[0].plot(timeScaleToHrs*Btime,Bp,linestyle='-',color='b',linewidth=2,label="pressure")

      # Plot domain plastic strain
      Bepxx = jobObj.fields["Bepxx"].getData()
      Bepyy = jobObj.fields["Bepyy"].getData()
      Bepzz = jobObj.fields["Bepzz"].getData()
      Bevp = (Bepxx+Bepyy+Bepzz)

      Bphi_e = [ A*np.exp(- (p**equilibriumPorosityPressureExponent)/B) for p in Bp ]
      Bphi_p = [ np.exp(-e) * ( phi_i+np.exp(e) -1 ) for e in Bevp ]
      

      axes[3].plot(Btime,Bphi_e,linestyle='-',color='r',linewidth=3,alpha=0.5,label=filename+"_phi_e")
      axes[3].plot(Btime,Bphi_p,linestyle='-',color='b',linewidth=3,alpha=0.5,label=filename+"_phi_p")

      #axes[2].plot(Btime,-Bepxx,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepxx")
      #axes[2].plot(Btime,-Bepyy,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepyy")
      #axes[2].plot(Btime,-Bepzz,linestyle=':',color=colors[i],linewidth=1,label=filename+"_Bepzz")
      #axes[2].plot(Btime,-Bevp,linestyle='-',color='g',linewidth=3,alpha=0.5,label=filename+"_Bevp")
      #axes[2].plot(Btime,-Bev,linestyle='-',color='b',linewidth=3,alpha=0.5,label=filename+"_Bev")

    if plotReactions:
      engineeringStress = False
      jobObj.compute_domain_stress(engineeringStress)
      time = jobObj.fields["Time"].getData()
      
      # pot domain deformation
      F00 = jobObj.fields["F00"].getData()
      F11 = jobObj.fields["F11"].getData()
      F22 = jobObj.fields["F22"].getData()     
      
      ev = np.log(F00*F11*F22)
      ev = F00*F11*F22 - 1.
      axes[1].plot(time,F00,linestyle='-',color=colors[i],linewidth=1,label=filename+"_F00")
      axes[1].plot(time,F11,linestyle='--',color=colors[i],linewidth=1,label=filename+"_F11")
      axes[1].plot(time,F22,linestyle=':',color=colors[i],linewidth=1,label=filename+"_F22")
      

      #axes[2].plot( [0,.002,.0045,.007,.011,.0137,0.022],[0.225,0.224,] )

      
      # plot domain stress computed from reactions:
      Rsxx = jobObj.fields["Rsxx"].getData()
      Rsyy = jobObj.fields["Rsyy"].getData()
      Rszz = jobObj.fields["Rszz"].getData()
      # axes[0].plot(timeScaleToHrs*time,-Rsxx,linestyle='-',color='r',linewidth=1,label=filename+"_Rsxx")
      # axes[0].plot(timeScaleToHrs*time,-Rsyy,linestyle='-',color='g',linewidth=1,label=filename+"_Rsyy")
      # axes[0].plot(timeScaleToHrs*time,-Rszz,linestyle='-',color='b',linewidth=1,label=filename+"_Rszz")
      
      timeDataNorm=[time[-1] * i/timeData[-1] for i in timeData];
      
      
      strainDataShifted=[ 0.003 + i - strainData[500] for i in strainData];    
      
      # axes[2].plot(timeDataNorm,strainDataShifted,linestyle='-',color='b',linewidth=2,alpha=0.7,label="Data_ev")
      axes[2].plot(timeScaleToHrs*np.array(timeDataNorm[500:]),strainDataShifted[500:],linestyle='-',color='k',linewidth=2,alpha=1.0,label="data")
      axes[2].plot(timeScaleToHrs*np.array(time),-ev,linestyle='--',color='r',linewidth=2,label="model")
      
      # Rp = (-1./3.)*(Rsxx+Rsyy+Rszz)
      # p = interpolate.interp1d(list(time),list(Rp))
      # p_Dinterpolate = p(timeDataNorm)

      # elasticStrainData = [ i / 2. for i in p_Dinterpolate ]  # assumes bulk = b0 = 2GPa
      # plasticStrainData = np.array(strainDataShifted) - np.array(elasticStrainData)

      # phiDataShifted = [ np.exp(e) * ( phi_i+np.exp(-e) -1 ) for e in plasticStrainData ]
      # axes[3].plot(timeDataNorm,phiDataShifted,linestyle='-',color='g',linewidth=3,alpha=0.5,label=filename+"_phi_D")
            
     
    read_files.append(file)

  except BaseException:
    logging.exception("Could not read file: " + file)
    print()
    errored_files.append(file)

# Format plots:
plt.figure(1)
axes[0].set_ylabel(r'pressure (GPa)')
axes[0].tick_params(axis='both', which='major')
axes[0].tick_params(axis='both', which='minor')
axes[0].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[0].set_xlim(1,90)
axes[0].grid()

axes[1].set_ylabel(r'F')
axes[1].tick_params(axis='both', which='major')
axes[1].tick_params(axis='both', which='minor')
axes[1].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[1].grid()

axes[2].set_xlabel('time (h)')
axes[2].set_ylabel(r'vol. strain')
axes[2].tick_params(axis='both', which='major')
axes[2].tick_params(axis='both', which='minor')
axes[2].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[2].set_xlim(1,90)
axes[2].grid()

#axes[3].set_xlabel('Time (μs)')
axes[3].set_ylabel(r'poro')
axes[3].legend(bbox_to_anchor=(1.22,1), loc="upper left",fontsize=str(legend_font_size))
axes[3].grid()

fig.tight_layout()
fig.savefig("geomechanicsHydrostaticCreep.png", bbox_inches="tight")

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
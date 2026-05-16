# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Geometry object functions for the particle file writer.
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
import importlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from os.path import exists
from cycler import cycler
import csv

fig, axes = plt.subplots(2, 2, figsize=(12,12) )

def lighten_color(color, amount=0.5):
  """
  Lightens the given color by multiplying (1-luminosity) by the given amount.
  Input can be matplotlib color string, hex string, or RGB tuple.
  Examples:
  >> lighten_color('g', 0.3)
  >> lighten_color('#F034A3', 0.6)
  >> lighten_color((.3,.55,.1), 0.5)
  """
  import matplotlib.colors as mc
  import colorsys
  try:
    c = mc.cnames[color]
  except:
    c = color
  c = colorsys.rgb_to_hls(*mc.to_rgb(c))
  return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

def filter(x):
  window_size=2
  x = np.array(x)
  N = len(x)
  x_out = np.copy(x)
  for n in range(window_size+1, N-window_size):
    x_out[n] = np.median(x[(n-window_size):(n+window_size+1)])
  return x_out 

# Stress strain data for which fitting parameters in this file were obtained in a stand-alone Mathematic notebook.
dataTension_strain=np.array([0.,0.0152041,0.0257091,0.0399285,0.0613519,0.0899947,0.125734,0.190025,0.254301,0.404224,0.586254,0.768264,0.960942,1.06796,1.06986])
dataTension_stress=np.array([0.1,17.1782,25.7759,29.3165,29.6106,25.7825,24.5412,24.3057,24.7979,28.7834,35.8106,46.3524,65.5892,83.2086,2.09566])

dataCompression_strain=np.array([0,0.00456397,0.012768,0.0249448,0.0471481,0.0683025,0.0874232,0.126676,0.165935,0.21224,0.326003,0.386414,0.44078,0.500182,0.500155,0.479909,0.462677,0.450462,0.446298])
dataCompression_stress=np.array([0.1,3.98692,11.183,21.6455,31.7168,34.7841,33.196,31.3032,30.5627,29.8357,30.8236,32.6245,33.3545,34.8972,28.9904,13.4961,5.99795,2.37367,0.9284])

axes[0][0].plot(-1.0*dataCompression_strain,-1.0*dataCompression_stress,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Data - Compression')
axes[0][0].plot(dataTension_strain,dataTension_stress,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Data - Tension')
        
# parent directory of file pathL:
runLocation='/p/lustre1/homel1/geosxRuns/'

files=[
'polymerCompression',
'polymerTension',
]  

labels=[
'Baseline compression response, no Temp Dependence',
'Baseline tension response, no Temp Dependence'
]  

styles=[
'-',
'-',
'--',
'-',
'--',
'-',
'--',
'-'
]

#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]

good=[]
for i,file in enumerate(files):
  print(runLocation+file+'/reactionHistory.csv')
  if ( exists(runLocation+file+'/reactionHistory.csv') ):
    good.append(i)

print("good = ",good)

for i in good:
  file = files[i]
  print(file)
  sys.path.append(runLocation+file+"/")
  print('Importing job')
  job = importlib.import_module('pfw_input_'+file)
  print('Job Imported.')
  pfw = job.pfw
  label = labels[i]
  reactionFile=runLocation+file+'/reactionHistory.csv'
    

  # Reaction from simulation results
  data = np.genfromtxt(reactionFile, delimiter=',')
  time = data[:,0]
  
  # Box values from simulation results.
  boxFile=runLocation+file+'/boxAverageHistory.csv'
  boxdata = np.genfromtxt(boxFile, delimiter=',')
  box_time = boxdata[1:,0]
  
  # for incomplete writes, lengths may be different for box and reaction files.  select shorter and truncate:
  l0=min(len(box_time),len(time))
  time = data[1:l0,0]
  box_time = boxdata[1:l0,0]

  # Cut reaction arrays to length
  F00 = data[1:l0,1]
  F11 = data[1:l0,2]
  F22 = data[1:l0,3]
  length_x = data[1:l0,4]
  length_y = data[1:l0,5]
  length_z = data[1:l0,6]
  Rxm = data[1:l0,7]
  Rxp = data[1:l0,8]
  Rym = data[1:l0,9]
  Ryp = data[1:l0,10]
  Rzm = data[1:l0,11]
  Rzp = data[1:l0,12]
  
  # Cut box arrays to length:
  box_sxx = boxdata[1:l0,1]
  box_syy = boxdata[1:l0,2]
  box_szz = boxdata[1:l0,3]
  box_rho = boxdata[1:l0,7]
  box_damage = boxdata[1:l0,8]
  box_e = boxdata[1:l0,9]
  box_ke = boxdata[1:l0,10]
  box_epxx = np.array(boxdata[1:l0,11])
  box_epyy = np.array(boxdata[1:l0,12])
  box_epzz = np.array(boxdata[1:l0,13])
  box_epyx = np.array(boxdata[1:l0,14])
  box_epxz = np.array(boxdata[1:l0,15])
  box_epxy = np.array(boxdata[1:l0,16])
  box_volume = boxdata[1:l0,17]
  box_temperature = boxdata[1:l0,18]
  box_F00 = boxdata[1:l0,19]
  box_F11 = boxdata[1:l0,20]
  box_F22 = boxdata[1:l0,21]
  
  # this hides all non-monotic time entries, so restart data files are cleaned:
  maxt = 0.0
  mask = np.ones(len(time), dtype=bool)
  avedt = time[-1]/len(time)
  for ii,t in enumerate(time):
    if ( t<=maxt ):
      mask[ii] = False
    else:
      maxt = t

  time = time[mask,...]
  F00 = F00[mask,...]
  F11 = F11[mask,...]
  F22 = F22[mask,...]
  Rxm = Rxm[mask,...]
  Rxp = Rxp[mask,...]
  Rym = Rym[mask,...]
  Ryp = Ryp[mask,...]
  Rzm = Rzm[mask,...]
  Rzp = Rzp[mask,...]
  length_x = length_x[mask,...]
  length_y = length_y[mask,...]
  length_z = length_z[mask,...]

  box_time = box_time[mask,...]
  box_sxx = box_sxx[mask,...]
  box_syy = box_syy[mask,...]
  box_szz = box_szz[mask,...]
  box_rho = box_rho[mask,...]
  box_e = box_e[mask,...]
  box_epxx = box_epxx[mask,...] 
  box_epyy = box_epyy[mask,...]
  box_epzz = box_epzz[mask,...]
  box_epyx = box_epyx[mask,...]
  box_epxz = box_epxz[mask,...]
  box_epxy = box_epxy[mask,...]
  
  box_damage = box_damage[mask,...]
  box_volume = box_volume[mask,...] 
  box_temperature = box_temperature[mask,...] 
  box_F00 = box_F00[mask,...] 
  box_F11 = box_F11[mask,...] 
  box_F22 = box_F22[mask,...] 


  # reaction stress-strain
  area = np.array(box_volume / length_y) # Current cross section estimated from particle volume / sample length
  A0 = area[0]    # initial area from input file
  
  # Engineering stress/initial area
  syy0 = 0.5*(Ryp-Rym) / A0
  # True stress/current area
  syy = 0.5*(Ryp-Rym) / area
  # "true" logarithmic strain:
  eyy = np.log(F11)
   
  # plot reaction stress values"
  axes[0][0].plot(filter(eyy), -filter(-1000.*syy),linestyle='--',color=colors[i],linewidth=2,label=label+r' $(R_y/A)$')
  axes[0][0].plot(filter(eyy), -filter(-1000.*syy0),linestyle=':',color=colors[i],linewidth=2,label=label+r' $(R_y/A_0)$')

  # box stress-strain
  box_eyy=np.log(box_F11)
  # box_syy is continuum stress averaged over domain, we must correct for solid volume fraction:
  box_solidVolumeFraction = box_volume / (length_x * length_y * length_z)
  box_syyMaterial = box_syy / box_solidVolumeFraction
  axes[0][0].plot(filter(box_eyy),-filter(-1000.*box_syyMaterial),linestyle='-',color=colors[i],linewidth=2,label=label+' (box)')
  
  # Plot equivalent plastic strain:
  box_eqep = np.sqrt((4./9)*0.5*((box_epxx - box_epyy)**2 + (box_epyy - box_epzz)**2 + (box_epzz - box_epxx)**2 + 6*(box_epxy**2 + box_epxz**2 + box_epyx**2)))

  axes[0][1].plot(filter(box_eyy),filter(box_eqep),linestyle='--',color=colors[i],linewidth=2,label=label)
  
  # Plot Damage evolution: 
  axes[1][0].plot(filter(box_eyy),filter(box_damage),linestyle='--',color=colors[i],linewidth=2,label=label)
   
  # Plot temperature
  axes[1][1].plot(filter(box_eyy),filter(box_temperature),linestyle='--',color=colors[i],linewidth=2,label=label)


# Format stress-strain plot:
axes[0][0].set_ylim(-100.,100.)
axes[0][0].set_xlim(-0.5,1.25)
axes[0][0].grid()
axes[0][0].legend()
axes[0][0].set_xlabel(r'strain (mm)', fontsize=16)
axes[0][0].set_ylabel(r'Compressive stress (MPa)', fontsize=16)
axes[0][0].tick_params(axis='both', which='major', labelsize=12)

# Format Damage plot
axes[1][0].set_xlabel(r'strain', fontsize=16)
axes[1][0].set_ylabel(r'Average Damage', fontsize=16)

# Format Plastic strain plot
axes[0][1].set_xlabel(r'strains', fontsize=16)
axes[0][1].set_ylabel(r'Average Equiv. P. Strain', fontsize=16)

# Format temperature Plot
axes[1][1].set_xlabel(r'strains', fontsize=16)
axes[1][1].set_ylabel(r'Temperature (K)', fontsize=16)

fig.tight_layout()
fig.savefig("polymerCompressionTension.png", bbox_inches="tight")
      # -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Modified for CDTriaxial plotting
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt
import sys
import matplotlib.cm as cm
import importlib
import matplotlib.pyplot as plt
from cycler import cycler
import csv


# def lighten_color(color, amount=0.5):
#     """
#     Lightens the given color by multiplying (1-luminosity) by the given amount.
#     Input can be matplotlib color string, hex string, or RGB tuple.

#     Examples:
#     >> lighten_color('g', 0.3)
#     >> lighten_color('#F034A3', 0.6)
#     >> lighten_color((.3,.55,.1), 0.5)
#     """
#     import matplotlib.colors as mc
#     import colorsys
#     try:
#         c = mc.cnames[color]
#     except:
#         c = color
#     c = colorsys.rgb_to_hls(*mc.to_rgb(c))
#     return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])

# Sometimes the reaction data are noisy, this cleans them up.
# box average data shouldn't be noisy, if it is, there is something
# actually wrong with the results.
def filter(x):
	window_size=2
	x = np.array(x)
	N = len(x)
	x_out = np.copy(x)
	for n in range(window_size+1, N-window_size):
		x_out[n] = np.median(x[(n-window_size):(n+window_size+1)])
	return x_out 

# parent directory of file path - updated for CDTriaxial:
runLocation='/data1/sghosh29/Working_MPM_LLNL/testGEOS/CDtriaxial/'
file='CDtriaxial'  # Single file for CDTriaxial

label='CDTriaxial'


sampleHeight = 1. # mm
sampleWidth = 1.  # mm

style='-'

#initialize plots
color = cm.rainbow(0.5)  # Use a single color for CDTriaxial
fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2,figsize=(24,24))
fig.suptitle('CDTriaxial Compression Test - Ceramic Damage Model',fontsize='large')

# BOX AVERAGE HISTORY FILE -------------------------------------------------------------------------------
# plot box sum data and compute initial volume fraction
print(file)
sys.path.append(runLocation)
job = importlib.import_module('pfw_input_'+file)

# Get material properties from input file
density = job.density 
bulk = job.bulk
shear = job.shear
tensileStrength = job.tensileStrength
compressiveStrength = job.compressiveStrength
maximumStrength = job.maximumStrength
crackSpeed = job.crackSpeed
confiningPressure = job.confiningPressure
stopTime = job.stopTime
maxCompressiveStrain = job.maxCompressiveStrain

reactionFile=runLocation+'boxAverageHistory.csv'

data = np.genfromtxt(reactionFile, delimiter=',')
box_time = data[:,0]
box_avg_Sxx = data[:,1]
box_avg_Syy = data[:,2]
box_avg_Szz = data[:,3]
box_avg_exx = data[:,11]  # these are plastic strains
box_avg_eyy = data[:,12]
box_avg_ezz = data[:,13]

box_F00 = data[:,18]
box_F11 = data[:,19]
box_F22 = data[:,20]

damage = data[:,8]

maxt = 0.0
mask = np.ones(len(box_time), dtype=bool)
for ii,t in enumerate(box_time):
	if (t<=maxt):
		mask[ii] = False
	else:
		maxt = t
box_time = box_time[mask,...]
box_avg_Sxx = box_avg_Sxx[mask,...]
box_avg_Syy = box_avg_Syy[mask,...]
box_avg_Szz = box_avg_Szz[mask,...]
box_avg_exx = box_avg_exx[mask,...]
box_avg_eyy = box_avg_eyy[mask,...]
box_avg_ezz = box_avg_ezz[mask,...]

# Print header for clarity
print(f"{'Index':>5} {'Sxx':>12} {'Syy':>12} {'Szz':>12}")

# Loop over the first 10 values and print them in one row each
for i, (sxx, syy, szz) in enumerate(zip(box_avg_Sxx[:20], box_avg_Syy[:20], box_avg_Szz[:20])):
    print(f"{i:5d} {sxx:12.6f} {syy:12.6f} {szz:12.6f}")



box_F00 = box_F00[mask,...]
box_F11 = box_F11[mask,...]
box_F22 = box_F22[mask,...]

box_exx = np.array(box_F00) - 1.
box_eyy = np.array(box_F11) - 1.
box_ezz = np.array(box_F22) - 1.

box_avg_diff_stress = box_avg_Sxx - box_avg_Szz	
box_avg_vol_stress = (-1.0/3.0)*(box_avg_Sxx+box_avg_Syy+box_avg_Szz)
box_avg_vol_strain = (-1.0)*(box_avg_exx+box_avg_eyy+box_avg_ezz)
print('job progress: ',box_time[-1]/stopTime)

# REACTION HISTORY FILE -------------------------------------------------------------------------------
reactionFile=runLocation+'reactionHistory.csv'
# area (mm^2)
A0=sampleWidth**2

# time,F00,F11,F22,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
data = np.genfromtxt(reactionFile, delimiter=',')
time = data[1:,0]
F00 = data[1:,1]
F11 = data[1:,2]
F22 = data[1:,3]
lx  = data[1:,4]
ly  = data[1:,5]
lz  = data[1:,6]
#R is reaction if kvarious directions
# m = minus and p = +
Rxm = data[1:,7]
Rxp = data[1:,8]
Rym = data[1:,9]
Ryp = data[1:,10]
Rzm = data[1:,11]
Rzp = data[1:,12]

# this hides all non-monotic time entries, so restart data files are cleaned:
maxt = 0.0
mask = np.ones(len(time), dtype=bool)
for ii,t in enumerate(time):
	if (t<=maxt):
		mask[ii] = False
	else:
		maxt = t
time =filter(time[mask,...])
F00 = filter(F00[mask,...])
F11 = filter(F11[mask,...])
F22 = filter(F22[mask,...])
lx = filter(lx[mask,...])
ly = filter(ly[mask,...])
lz = filter(lz[mask,...])
Rxm = filter(Rxm[mask,...])
Rxp = filter(Rxp[mask,...])
Rym = filter(Rym[mask,...])
Ryp = filter(Ryp[mask,...])
Rzm = filter(Rzm[mask,...])
Rzp = filter(Rzp[mask,...])

# strain
exx=np.log(F00)
eyy=np.log(F11)
ezz=np.log(F22)
ev = (-1. )*(exx + eyy + ezz)

#stress
Ax=ly*lz
Ay=lx*lz
Az=lx*ly

#stress in x, y z directions
sxx=0.5*(Rxp-Rxm)/Ax
syy=0.5*(Ryp-Rym)/Ay
szz=0.5*(Rzp-Rzm)/Az

rxn_diff_stress=szz-sxx

#p= volumetric or hydrostatic stress
p=(-1.0/3.0)*(sxx+syy+szz)
rho=density/F00

# stress is scaled from simulation units (GPa) to plotting units (MPa)
ax1.plot(time,F00,linestyle='-',color=cm.gist_rainbow(0),linewidth=1)
ax1.plot(time,F11,linestyle='-',color=cm.gist_rainbow(0.25),linewidth=1)
ax1.plot(time,F22,linestyle='-',color=cm.gist_rainbow(0.75),linewidth=1)

# PLOTTING SIMULATED DATA -------------------------------------------------------------------------------

# plot model results
ax2.plot(box_avg_vol_stress, box_avg_diff_stress,linestyle=style,color=color,linewidth=2,label='CDTriaxial model')

#plot differential stress versus axial and radial strain
ax4.plot(-box_exx,box_avg_diff_stress,linestyle=':',color=color,linewidth=2,label=label+' model lat strain')
ax4.plot(-box_ezz,box_avg_diff_stress,linestyle='--',color=color,linewidth=2,label=label+' model ax strain')

#plot reaction stress and box average stress versus time 
ax5.plot(box_time,box_avg_vol_stress,linestyle=style,color=cm.gist_rainbow(0.2),linewidth=2,label='box averaged volumetric stress')

# Plot damage evolution
ax6.plot(box_time, damage, linestyle=style, color=color, linewidth=2, label='Damage evolution')

#plot Sxx, Syy, Szz vs time in ax3
ax3.plot(box_time, -box_avg_Sxx, linestyle='-', color='r', linewidth=2, label='Sxx')
ax3.plot(box_time, -box_avg_Syy, linestyle='--', color='g', linewidth=2, label='Syy')
ax3.plot(box_time, -box_avg_Szz, linestyle=':', color='b', linewidth=2, label='Szz')

# LABELING/FORMATTING PLOTS -------------------------------------------------------------------------------

ax1.set_xlabel(r'time ($\mu$s) ', fontsize=16)
ax1.set_ylabel('F00,F11,F22', fontsize=16)
# No legend for ax1 as per user request
ax1.grid()

ax2.set_xlabel(r'$pressure$ (GPa)', fontsize=16)
ax2.set_ylabel(r'$differential stress$ (GPa)', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)
handles, labels = ax2.get_legend_handles_labels()
if labels and any(label and not label.startswith('_') for label in labels):
    ax2.legend(handles, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize='medium')
ax2.grid()

ax3.set_xlabel('time ($\mu$s)', fontsize=16)
ax3.set_ylabel(r'$-\sigma_{xx,yy,zz}$ (GPa)', fontsize=16)
handles, labels = ax3.get_legend_handles_labels()
if labels and any(label and not label.startswith('_') for label in labels):
    ax3.legend(handles, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize='medium')
ax3.grid()

ax4.set_xlabel('$-\epsilon_{axial and lateral}$', fontsize=16)
ax4.set_ylabel(r'$-\sigma_{differential}$ (GPa)', fontsize=16)
ax4.set_ylim(0.,0.045)
handles, labels = ax4.get_legend_handles_labels()
if labels and any(label and not label.startswith('_') for label in labels):
    ax4.legend(handles, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize='medium')
ax4.grid()

ax5.set_xlabel('time ($\mu$s)', fontsize=16)
ax5.set_ylabel(r'$-\sigma_{volumetric}$ (GPa)', fontsize=16)
handles, labels = ax5.get_legend_handles_labels()
if labels and any(label and not label.startswith('_') for label in labels):
    ax5.legend(handles, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize='medium')
ax5.grid()

ax6.set_xlabel('time ($\mu$s)', fontsize=16)
ax6.set_ylabel('Damage', fontsize=16)
handles, labels = ax6.get_legend_handles_labels()
if labels and any(label and not label.startswith('_') for label in labels):
    ax6.legend(handles, labels, bbox_to_anchor=(1.04,1), loc="upper left", fontsize='medium')
ax6.grid()

fig.tight_layout()
#plt.show()
fig.savefig("CDTriaxial_results_updated.png", bbox_inches="tight")


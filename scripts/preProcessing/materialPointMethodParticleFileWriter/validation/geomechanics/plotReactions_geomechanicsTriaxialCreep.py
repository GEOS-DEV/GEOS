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
import math
from cycler import cycler
import csv

import re

#for interpolating the time and strain data
from scipy import interpolate

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


# parent directory of file pathL:
runLocation='/p/lustre1/homel1/geosxRuns/'
#runLocation='/p/lustre1/malenda1/geosRuns/'
#runLocation='/usr/workspace/homel1/geosRuns/'

# Material model parameters.  Will be read from input file for plotting.
density = 2.648
bulk = 36.3
shear = 26.0
sigmaYield = 2.27
sigmaFail = 0.449
sigmaMax = 5.0
crackSpeed = 1.8

youngsModulus = (9.*bulk*shear)/(3.*bulk+shear)
strainToFailure = sigmaFail/youngsModulus
stretch = np.exp(2.0*strainToFailure) # stretch to 50% past the strain to failure.
weibullSeed=1
weibullModulus=6.0
weibullVolume=1.0
flawSize=0.2

sampleHeight = 1. # mm
sampleWidth = 1.  # mm


#replace the file names below with the run you just completed
files=[
'geoTXCr_2024_TC02',
#'strain_control_geoCreepHYD'
]  


labels=[
'TC02',
#'strain_control_geoCreepHYD'
]

testnames = ['TC02']
length_of_test_days = 2

confiningPressure = 0.01  # confining pressure GPa


styles=[
'-',
'--',
'--',
'--',
'--',
'--',
'--',
'-.',
'-.',
'-.'
]

stopTime = 1000.0

stressTable=np.array([[2.487418596875372521e-05*stopTime, 0, 0, 0],
[1.012563303257513030e-02*stopTime, -2.108803608999999975e-03, -2.108803608999999975e-03, -2.108803608999999975e-03],
[2.022639187918150866e-02*stopTime, -4.593658936999999755e-03, -4.593658936999999755e-03, -4.593658936999999755e-03],
[3.032715072578788529e-02*stopTime, -7.453454796999999833e-03, -7.453454796999999833e-03, -7.453454796999999833e-03],
[5.052866841900063855e-02*stopTime, -7.630466507000000821e-03, -7.630466507000000821e-03, -7.630466507000000821e-03],
[6.062942726560701517e-02*stopTime, -confiningPressure, - confiningPressure, -1.070101939000000016e-02],
[6.062942726560701517e-02*stopTime, -confiningPressure, - confiningPressure, -1.070101939000000016e-02],
[7.073018611221337792e-02*stopTime, -confiningPressure, - confiningPressure, -1.509734128000000085e-02],
[9.093170380542614506e-02*stopTime, -confiningPressure, - confiningPressure, -1.547444327999999980e-02],
[1.010324626520325147e-01*stopTime, -confiningPressure, - confiningPressure, -1.795881333000000010e-02],
[1.010324626520325147e-01*stopTime, -confiningPressure, - confiningPressure, -1.795881333000000010e-02],
[1.111332214986388844e-01*stopTime, -confiningPressure, - confiningPressure, -2.113893686000000230e-02],
[1.313347391918516516e-01*stopTime, -confiningPressure, - confiningPressure, -2.104603639999999901e-02],
[1.414354980384580074e-01*stopTime, -confiningPressure, - confiningPressure, -2.485571896000000197e-02],
[1.414354980384580074e-01*stopTime, -confiningPressure, - confiningPressure, -2.485571896000000197e-02],
[stopTime, -confiningPressure, - confiningPressure, -2.615328664999999955e-02]])

################################# initialize plots ##################################################################################
#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
#fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(16,16))
fig, ((ax1, ax2,ax3, ax4)) = plt.subplots(4, 1,figsize=(16,16))

################################# get informatoin from input file                                       #################################
################################# the strain are plotted relative to the hydrostatic state at t=tStop/2 #################################
################################# (change this if stress table time changes)                            #################################
for i,file in enumerate(files):	 
	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)
 
	confiningPressure = job.confiningPressure

	stren = job.STREN
	peakI1 = job.PEAKI1
	ySlope = job.YSLOPE
	fSlope = job.FSLOPE
	fSlopeFailed = job.FSLOPEFAILED
	strainHardeningK = job.strainHardeningK
	strainHardeningn = job.strainHardeningn
	fractureEnergyReleaseRate = job.fractureEnergyReleaseRate
	fractureStress = job.fractureStress
	b0 = job.b0
	beta = job.BETA_nonassociativity
	g0 = job.g0
	g1 = job.g1
	g2 = job.g2
 
	c0 = job.creepc0
	c1 = job.creepc1
	c2 = job.creepc2
	D = job.creepD
	E = job.creepE
	F = job.creepF

	CR = job.CR
	p0 = job.p0
	a1 = stren
	a2 = (fSlope-ySlope )/(stren-ySlope*peakI1)
	a3 = (stren-ySlope*peakI1)*np.exp(-a2*peakI1)
	a4 = ySlope

	a1f = stren
	a2f = (fSlopeFailed-ySlope )/(stren)
	a3f = stren
	a4f = ySlope

	strenh = (stren+strainHardeningK)
	peakI1h = peakI1 + strainHardeningK/fSlope
	a1h = strenh
	a2h = (fSlope-ySlope )/(strenh-ySlope*peakI1h)
	a3h = (strenh-ySlope*peakI1h)*np.exp(-a2*peakI1h)
	a4h = ySlope
 
	stopTime=job.stopTime
	timeScale=job.timeScale*1e-6 # scales time to seconds

 
	print('g0='+str(g0), 'c0='+str(c0), 'c1='+str(c1), 'c2='+str(c2))
 

################################## determine original yield envelope values ################################# 
#	I1plot = np.linspace(p0,peakI1,100)
#	pPlot = np.array([ -I1/3. for I1 in I1plot ])
#	Ff_MPa = 1000.*np.array([Ff(I1,a1,a2,a3,a4) for I1 in I1plot ])
#	FfFc_MPa = Ff_MPa*np.array([Fc(I1,p0,peakI1,CR) for I1 in I1plot ])
# 
################################## determine damage (or failed) yield envelope values #################################  
#	I1Fplot = np.linspace(p0,0,100)
#	pFPlot = np.array([ -I1/3. for I1 in I1Fplot ])
#	FfFailed_MPa = 1000.*np.array([Ff(I1,a1f,a2f,a3f,a4f) for I1 in I1Fplot ])
#	FfFailedFc_MPa = FfFailed_MPa*np.array([Fc(I1,p0,0,CR) for I1 in I1Fplot ])
# 
################################## determine hardened (or failed) yield envelope values #################################  
#	I1Hplot = np.linspace(p0,peakI1h,100)
#	pHPlot = np.array([ -I1/3. for I1 in I1Hplot ])
#	FfHardened_MPa = 1000.*np.array([Ff(I1,a1h,a2h,a3h,a4h) for I1 in I1Hplot ])
#	FfHardenedFc_MPa = FfHardened_MPa*np.array([Fc(I1,p0,peakI1h,CR) for I1 in I1Hplot ])	
# 
 

 
################################# pull box data  ################################################################## 
 
	boxDataFile=runLocation+file+'/boxAverageHistory.csv'
	# Time, Sxx, Syy, Szz, Syz, Sxz, Sxy, Density, Damage, Internal Energy, Kinetic Energy, epxx, epyy, epzz, epyz, epxz, epxy, volume
	data = np.genfromtxt(boxDataFile, delimiter=',')
	if(data.ndim < 2):
		print('bad data file')
	else:
		box_time = data[1:,0]
		box_sxx = data[1:,1]
		box_syy = data[1:,2]
		box_szz = data[1:,3]
		box_sxy = data[1:,4]
		box_syz = data[1:,5]
		box_sxz = data[1:,6]
		box_rho = data[1:,7]
		box_damage = data[1:,8]
		box_e = data[1:,9]
		box_ke = data[1:,10]
		box_epxx = data[1:,11]
		box_epyy = data[1:,12]
		box_epzz = data[1:,13]
		# box_epyz = data[1:,14]
		# box_epxz = data[1:,15]
		# box_epxy = data[1:,16]
		box_volFrac = data[1:,17]
		box_F00 = data[1:,18]
		box_F11 = data[1:,19]
		box_F22 = data[1:,20]
  
  
		#print('box data',box_time, box_sxx)# this hides all non-monotic time entries, so restart data files are cleaned:
		maxt = 0.0
		t = 0.0
		mask = np.ones(len(box_time), dtype=bool)
		for ii,t in enumerate(box_time):
			if (t<=maxt):
				mask[ii] = False
			else:
				maxt = t
		box_time = box_time[mask,...]
		box_sxx = box_sxx[mask,...]
		box_syy = box_syy[mask,...]
		box_szz = box_szz[mask,...]
		box_sdiff = box_szz - box_sxx
		box_svol = (-1.0/3.0)*(box_sxx+box_syy+box_szz)
		box_sxy = box_sxy[mask,...]
		box_syz = box_syz[mask,...]
		box_sxz = box_sxz[mask,...]
		box_rho = box_rho[mask,...]
		box_epxx = box_epxx[mask,...]  
		box_epyy = box_epyy[mask,...]  
		box_epzz = box_epzz[mask,...]
		# box_epyz = box_epyz[mask,...]  
		# box_epxz = box_epxz[mask,...]  
		# box_epxy = box_epxy[mask,...]
	
		box_e = box_e[mask,...]
		box_damage = box_damage[mask,...]
		box_volFrac = box_volFrac[mask,...]
		box_F00 = box_F00[mask,...]  
		box_F11 = box_F11[mask,...]  
		box_F22 = box_F22[mask,...]
  
		box_exx = (box_F00) - 1
		box_eyy = (box_F11) - 1
		box_ezz = (box_F22) - 1
		box_evol = box_exx + box_eyy+ box_ezz
 
 
		box_p = (-1.0/3.0)*(box_sxx+box_syy+box_szz)
		box_vm = np.sqrt(0.5*( np.square(box_sxx-box_syy)+np.square(box_syy-box_szz)+np.square(box_szz-box_sxx) ) )
		box_ev = -np.log(box_rho[0] / np.array(box_rho))
		box_evp = box_epxx + box_epyy + box_epzz
		#box_vm_shear = np.sqrt(3)*box_sxz
		equivalent_shear_strain = (2/3)*np.sqrt((box_exx**2)-(box_exx*(box_eyy+box_ezz))+(box_eyy**2)-(box_eyy*box_ezz)+(box_ezz**2))
  
		equilibriumShearStrainConstant = c0
		equilibriumShearStrainExponent = c1
		ShearStrainRateConstant = c2
  
		#find J2
		Tau = []
		J2 = []
		equilibriumVMShearStrain = []
		for i in range(0, len(box_szz)):
			Tau.append((box_sxx[i]-box_szz[i])/(np.sqrt(3)))
			#J2.append(Tau[i]**2)
			#in the hpp file, plastc VM shear strain is: real64 plasticVMshearStrain = sqrt (2./3.) * ep_rootJ2_old;
			equilibriumVMShearStrain.append(equilibriumShearStrainConstant * (Tau[i] **( equilibriumShearStrainExponent)))
   
   


################################# GET EXPERIMENTAL	STRESS DATA ##################################################################
colorvec2=[0.3, 0.5, 0.8]	
colorval2=0
#
#
#
Tau_experiments = []
for test in testnames:
	# Open and read the CSV file
	with open('make_stress_table_TC02.csv', mode='r',encoding='utf-8-sig') as csv_file:  # Replace 'data.csv' with your filename
		reader = csv.reader(csv_file)
		data = list(reader)
		# will need to check individual excel files to see what these indices should be. 
		#short_data = data[0:(2900-900)]
	
		# Extract columns as lists
		experimental_stress_time = [float(row[0]) for row in data]  # Column 0
		experimental_stress_lat = [float(row[1])/1000 for row in data]  # Column 0
		experimental_stress_ax = [float(row[2])/1000 for row in data]  # Column 1
		experimental_stress_time = [(x*stopTime/experimental_stress_time[-1])for x in experimental_stress_time]


for i in range(0, len(experimental_stress_ax)):
    Tau_experiments.append((experimental_stress_ax[i] - experimental_stress_lat[i])/(np.sqrt(3.)))


################################## GET EXPERIMENTAL	STRAIN DATA ##################################################################
#
for test in testnames:
	# Open and read the CSV file
	with open('TC02_vol_ax_lat_strain.csv', mode='r',encoding='utf-8-sig') as csv_file:  # Replace 'data.csv' with your filename
		reader = csv.reader(csv_file)
		data = list(reader)
		# will need to check individual excel files to see what these indices should be. 
		
	
		# Extract columns as lists
		experimental_strain_time = [float(row[0]) for row in data]  # Column 2 (fewer rows)

		experimental_strain_vol = [float(row[1]) for row in data ]  # Column 2 (fewer rows)
		experimental_strain_ax = [float(row[2]) for row in data ]  # Column 2 (fewer rows)
		experimental_strain_lat = [float(row[3]) for row in data ]  # Column 2 (fewer rows)

		experimental_strain_time = [(x*stopTime/experimental_strain_time[-1])for x in experimental_strain_time]




################################# plot simulation and model data ################################################################## 
	#F00, F11, F22 versus time

#MH comments:
# In evaluating the fit, I'd probably create one plot that is lateral, axial, and volumetric strain (total not plastic) vs time 
# to compare directly with the experimental data.  

for i,file in enumerate(files):	 
	ax1.plot(timeScale*np.array(experimental_strain_time),experimental_strain_vol,linestyle='--',color=cm.Greys(0.25),linewidth=4,label='experimental volumetric strain')
	ax1.plot(timeScale*np.array(experimental_strain_time),experimental_strain_ax,linestyle=styles[i],color=cm.Greys(0.4),linewidth=4,label='experimental axial strain')
	ax1.plot(timeScale*np.array(experimental_strain_time),experimental_strain_lat,linestyle=styles[i],color=cm.Greys(0.9),linewidth=4,label='experimental lateral strain')

	ax1.plot(timeScale*np.array(box_time),-box_evol,linestyle='--',color=cm.Reds(0.25),linewidth=2,label=testnames[i] +' box strain volumetric')
	ax1.plot(timeScale*np.array(box_time),-box_exx,linestyle=styles[i],color=cm.Reds(0.25),linewidth=2,label=testnames[i] +' box strain xx ')
	ax1.plot(timeScale*np.array(box_time),-box_eyy,linestyle=styles[i],color=cm.Reds(0.9),linewidth=2,label=testnames[i] +' box strain yy')
	ax1.plot(timeScale*np.array(box_time),-box_ezz,linestyle=styles[i],color=cm.Reds(0.4),linewidth=2,label=testnames[i] +' box strain zz')



	#I'd have another plot that is stress components vs. time, so you can make sure the stress-control is working, 
	# 2) plot of x, y and z stresses versus time
	ax2.plot(timeScale*np.array(experimental_strain_time),experimental_stress_lat,linestyle=styles[i],color='k',linewidth=3,label='experimental lateral stress')
	ax2.plot(timeScale*np.array(experimental_strain_time),experimental_stress_ax,linestyle=styles[i],color=cm.Greys(0.5),linewidth=3,label='experimental axial stress')
	ax2.plot(timeScale*np.array(box_time),-box_sxx,linestyle=styles[i],color=cm.Greens(0.25),linewidth=3,label='-box x stress')
	ax2.plot(timeScale*np.array(box_time),-box_syy,linestyle=styles[i+1],color=cm.Greens(0.8),linewidth=3,label='-box y stress')
	ax2.plot(timeScale*np.array(box_time),-box_szz,linestyle=styles[i],color=cm.Greens(0.4),linewidth=3,label='-box z stress')
 
	print('Run completion = ',box_time[-1]/stopTime)


	#and I'd have another showing von Mises shear stress vs time,
	# 3) plot of vonmisses stresses versus time
	ax3.plot(timeScale*np.array(box_time),Tau, linestyle=styles[i],linewidth = 3, color=cm.Blues(0.5),label=' Tau')
	ax3.plot(timeScale*np.array(experimental_strain_time),Tau_experiments, linestyle=styles[i],linewidth = 2, color=cm.Greys(0.5),label=' Tau from experiments')
	

	#and one showing equivalent shear strain vs time (along with the data), which will be easier to interpret than looking at just the components
	# 4) plot of equivalent shear strain versus time
	ax4.plot(timeScale*np.array(experimental_strain_time),0.816497*(np.array(experimental_strain_ax)-np.array(experimental_strain_lat)),linestyle='-',color=cm.Greys(0.9),linewidth=6,label='Data: strain difference')
	ax4.plot(timeScale*np.array(box_time),-0.816497*(np.array(box_ezz)-np.array(box_exx)),linestyle='--',linewidth = 4, color='r',label='Model: strain difference')
	ax4.plot(timeScale*np.array(box_time),-0.816497*(np.array(box_epzz)-np.array(box_epxx)),linestyle='-',linewidth = 4, color='m',label='Model: plastic strain difference')
	ax4.plot(timeScale*np.array(box_time),equilibriumVMShearStrain,linestyle='--',linewidth = 4, color='b',label='Model: Equilibrium plastic shear strain')
	








################################# Labeling and formatting plots ##################################################################


ax1.set_xlabel(r'time (s)', fontsize=16)
ax1.set_ylabel(r'Strain (-)', fontsize=16)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax1.grid()
ax1.set_xlim(0, 200*timeScale)
ax1.set_ylim(0, 0.02)


ax2.set_xlabel(r'time (s)', fontsize=16)
ax2.set_ylabel(r'stress components (GPa)', fontsize=16)
ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax2.grid()
ax2.set_xlim(0, 200*timeScale)
ax3.set_ylim(0, 0.03)



ax3.set_xlabel(r'time (s)', fontsize=16)
ax3.set_ylabel(r'von Mises shear stress (GPa)', fontsize=16)
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax3.grid()
ax3.set_xlim(0, 200*timeScale)
ax3.set_ylim(0, 0.012)


ax4.set_xlabel(r'time (s)', fontsize=16)
ax4.set_ylabel(r'equivalent shear strain', fontsize=16)
ax4.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax4.grid()
ax4.set_xlim(0, 200*timeScale)
ax4.set_ylim(0, 0.035)




# SAVE PLOTS -------------------------------------------------------------------------------

fig.suptitle('TC02: c0'+str(c0)+ ' c1='+str(c1)+' c2=' +str(round(c2,3))+'D=' +str(D)+' E=' +str(E)+' F=' +str(F)+' g0=' +str(g0)+' b0=' +str(round(b0,3)),fontsize=20)

#fig.suptitle('ISR03 CR'+str(CR)+ ' k='+str(strainHardeningK)+' n=' +str(strainHardeningn) +'PEAKI1='+str(peakI1)+'STREN='+str(stren)+' \n Fslope='+str(fSlope)+ ' Yslope='+str(ySlope)+' p0='+str(p0)+' b0='+str(b0)+'beta='+str(beta)+'FERR='+str(fractureEnergyReleaseRate)+'frcstss='+str(round(fractureStress,3))+'g0='+str(round(g0,3))+'g1='+str(round(g1,3))+'g2='+str(round(g2,3)),fontsize=20)

fig.tight_layout()
#figname = "OTXCrD: c0:"+str(round(c0,4))+"CR:"+str(CR)+"_beta:"+str(beta)+"fsp:"+str(fSlope)+"ferr:"+str(fractureEnergyReleaseRate)+"hardK:"+str(strainHardeningK)+"hardn:"+str(strainHardeningn)+"stren:"+str(stren)+"g0="+str(round(g0,3))+"g1="+str(round(g1,3))+"g2="+str(round(g2,3))+'new9'
#figname = 'TC02: c0'+str(c0)+ ' c1='+str(c1)+' c2=' +str((round(c2,3)))+ 'D=' +str(D)+' E=' +str(E)+' F=' +str(F)+' g0=' +str((round(g0,3)))+' b0=' +str((round(b0,3)))

figname = 'geomechanicsTriaxialCreep_TC02'

#figname = "ISR03_V1_fracsts:"+str(round(fractureStress,4))+"CR:"+str(CR)+"_beta:"+str(beta)+"fsp:"+str(fSlope)+"ferr:"+str(fractureEnergyReleaseRate)+"hardK:"+str(strainHardeningK)+"hardn:"+str(strainHardeningn)+"stren:"+str(stren)+"g0="+str(round(g0,3))+"g1="+str(round(g1,3))+"g2="+str(round(g2,3))+'new9'
fig.savefig(figname+".png", bbox_inches="tight")
print(figname)
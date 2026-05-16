# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Geometry object functions for the particle file writer.
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt

import matplotlib.cm as cm
import sys

import importlib
import matplotlib.pyplot as plt
from cycler import cycler

import csv

#################################   Housekeeping      ################################################################
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


# parent directory of file pathL:
runLocation='/p/lustre1/homel1/geosxRuns/'
files=[
'geomechanicsTriaxialCompression_ISR03',
]  

testnames=[
'ISR03']

styles=[
'-',
'--',
'-.']

colorvec1=[0.3, 0.55, 0.8]	
colorval1=0

widths=[3.5,2.5,1.5]


################################# universal functions ##################################################################

def Ff(I1,a1,a2,a3,a4):
	return a1 - a3*np.exp(a2*I1)-a4*I1

def Fc(I1,p0,peakI1,CR):
	kappa = peakI1 + CR*(p0-peakI1)
	if (I1 >= kappa):
		return 1
	elif (I1 <= p0):
		return 0
	else:
		return np.sqrt(1 - ((-I1+kappa)/(-p0+kappa) )**2.)
	

################################# failure data from excel sheets  ##################################################################


#these next values are form the plots of the diff stress versus lateral and axial strain in the raw data excel sheets:
#/Users/malenda1/Desktop/creep-compaction/organized by test type/BTXCoD (ISR030506)/high temp
#lat_diff_stress_experiments = [32.39, 31.071, 40.8]
lat_diff_stress_experiments = [18., 18., 18.]
lat_diff_stress_experiments = [11., 18., 18.]

#lat_diff_stress_experiments = [item/1000 for item in lat_diff_stress_experiments]
#ax_diff_stress_experiments = [32.39, 31.071, 36.92]
ax_diff_stress_experiments = [14.44, 5.73, 14.14]
ax_diff_stress_experiments = [6.5, 5.73, 14.14]
#ax_diff_stress_experiments = [item/1000 for item in ax_diff_stress_experiments]


################################# initial material properties##################################################################

# input material properties
bulk = 36.6
shear = 26.0
C11 = bulk + (4.0/3.0)*shear
EE = 9.0*bulk*shear/(3.0*bulk + shear)
C12 = C11 - 2.0*shear


################################# initial volume fraction (assumes one material with density rho02)#################################
VF0=np.zeros(len(files))
rho0s=2.648

################################# initialize plots ##################################################################################
#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2,figsize=(16,16))


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



	print('stren='+str(stren)+'; peakI1='+str(peakI1)+'; ySlope='+str(ySlope)+'; fSlope='+str(fSlope)+'; fSlopeFailed='+str(fSlopeFailed)+'; Gf='+str(fractureEnergyReleaseRate)+'; sigmaF='+str(fractureStress)+'; g1='+str(g1)+'; g2='+str(g2))

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


 
	print('a1='+str(a1)+'; a2='+ str(a2) +'; a3='+ str(a3) +'; a4=' +str(a4))
 


################################# determine original yield envelope values ################################# 
	I1plot = np.linspace(p0,peakI1,100)
	pPlot = np.array([ -I1/3. for I1 in I1plot ])
	Ff_MPa = 1000.*np.array([Ff(I1,a1,a2,a3,a4) for I1 in I1plot ])
	FfFc_MPa = Ff_MPa*np.array([Fc(I1,p0,peakI1,CR) for I1 in I1plot ])
 
################################# determine damage (or failed) yield envelope values #################################  
	I1Fplot = np.linspace(p0,0,100)
	pFPlot = np.array([ -I1/3. for I1 in I1Fplot ])
	FfFailed_MPa = 1000.*np.array([Ff(I1,a1f,a2f,a3f,a4f) for I1 in I1Fplot ])
	FfFailedFc_MPa = FfFailed_MPa*np.array([Fc(I1,p0,0,CR) for I1 in I1Fplot ])
 
################################# determine hardened (or failed) yield envelope values #################################  
	I1Hplot = np.linspace(p0,peakI1h,100)
	pHPlot = np.array([ -I1/3. for I1 in I1Hplot ])
	FfHardened_MPa = 1000.*np.array([Ff(I1,a1h,a2h,a3h,a4h) for I1 in I1Hplot ])
	FfHardenedFc_MPa = FfHardened_MPa*np.array([Fc(I1,p0,peakI1h,CR) for I1 in I1Hplot ])	
 
################################# pull reaction data ##################################################################################
	rho0=rho0s*VF0[i]
	print(file)
	reactionFile=runLocation+file+'/reactionHistory.csv'
	# time,F00,F11,F22,lx,ly,lz,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
	# # time, F00, F11, F22, length_x, length_y, length_z, Rx-, Rx+, Ry-, Ry+, Rz-, Rz+, L00, L11, L22
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[:,0]
	F00 = data[:,1]
	F11 = data[:,2]
	F22 = data[:,3]
	lx  = data[:,4]
	ly  = data[:,5]
	lz  = data[:,6]
	Rxm = data[:,7]
	Rxp = data[:,8]
	Rym = data[:,9]
	Ryp = data[:,10]
	Rzm = data[:,11]
	Rzp = data[:,12]

	# this hides all non-monotic time entries, so restart data files are cleaned:
	# we then apply the filter to remove the spikes from the data do to instability at
	# driven boudnary.
	
	maxt = 0.0
	mask = np.ones(len(time), dtype=bool)
	for ii,t in enumerate(time):
		if (t<=maxt):
			mask[ii] = False
		else:
			maxt = t
	time = filter(time[mask,...])
	F00 = filter(F00[mask,...])
	F11 = filter(F11[mask,...])
	F22 =filter( F22[mask,...])
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
 
	ev = -np.log(F00*F11*F22)
 
 
	
	#stress
	Ax=ly*lz
	Ay=lx*lz
	Az=lx*ly

	sxx=0.5*(Rxp-Rxm)/Ax
	syy=0.5*(Ryp-Rym)/Ay
	szz=0.5*(Rzp-Rzm)/Az
 
	r_sdiff = szz - sxx

	p=(-1.0/3.0)*(sxx+syy+szz)
	rho=rho0/F00

	vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )

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
		# box_sxy = data[1:,4]
		# box_syz = data[1:,5]
		# box_sxz = data[1:,6]
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
		# box_sxy = box_sxy[mask,...]
		# box_syz = box_syz[mask,...]
		# box_sxz = box_sxz[mask,...]
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
  
		box_exx = np.log(box_F00)  
		box_eyy = np.log(box_F11)
		box_ezz = np.log(box_F22)
		box_evol = box_exx + box_eyy+ box_ezz
 
		box_p = (-1.0/3.0)*(box_sxx+box_syy+box_szz)
		box_vm = np.sqrt(0.5*( np.square(box_sxx-box_syy)+np.square(box_syy-box_szz)+np.square(box_szz-box_sxx) ) )
		box_ev = -np.log(box_rho[0] / np.array(box_rho))
		box_evp = box_epxx + box_epyy + box_epzz




################################# plot data ################################################################## 
	#F00, F11, F22 versus time
	ax1.plot(time,F00,linestyle=styles[i],color=cm.gist_rainbow(0),linewidth=1,label=testnames[i] +' reaction F00')
	ax1.plot(time,F11,linestyle=styles[i],color=cm.gist_rainbow(0.25),linewidth=1,label=testnames[i] +' reaction F11')
	ax1.plot(time,F22,linestyle=styles[i],color=cm.gist_rainbow(0.75),linewidth=1,label=testnames[i] +' reaction F22')
	ax1.plot(box_time, box_damage,linestyle='-',color='r',linewidth=1,label=testnames[i] +' damage')
 

	#ax2.plot(1000.*confiningPressure,(np.sqrt(1/3)*lat_diff_stress_experiments[i]),'o',color=cm.Greys(colorvec1[colorval1]),mec='k',markersize = 15, label=testnames[i]+' lateral experimental value')
	#ax2.plot(1000.*confiningPressure,(np.sqrt(1/3)*ax_diff_stress_experiments[i]),'^',color=cm.Greens(colorvec1[colorval1]),mec='k',markersize =10, label=testnames[i]+' axial experimental value')	
	ax2.plot(box_p,box_vm,linestyle=styles[i],color=cm.Oranges(colorvec1[colorval1]),linewidth=widths[i],label=testnames[i]+'  box vm')
	
	ax2.plot( 1000.*pPlot,FfFc_MPa,linestyle=styles[i],color=cm.Greens(colorvec1[colorval1]),linewidth=widths[i],label=testnames[i]+' Ff*Fc (non-hardening)')
	ax2.plot( 1000.*pFPlot,FfFailedFc_MPa,linestyle='--',color='r',linewidth=3,label='Fc*Ff_failed')	
	ax2.plot( 1000.*pHPlot,FfHardenedFc_MPa,linestyle=styles[i],color=cm.Blues(colorvec1[colorval1]),linewidth=widths[i],label=testnames[i]+' Ff*Fc (hardening)')	

	ax2.plot(1000.*box_p,np.sqrt(1/3)*1000.*box_vm,linestyle=styles[i],color=cm.Oranges(colorvec1[colorval1]),linewidth=widths[i],label=testnames[i]+'  box vm')


	#directly from Mike's input file
	ax3.plot(box_time,box_sxx,linestyle=styles[i],color=cm.Reds(0.1),linewidth=2,label='-box sxx')
	ax3.plot(box_time,box_syy,linestyle=styles[i],color=cm.Reds(0.4),linewidth=2,label='-box syy')
	ax3.plot(box_time,box_szz,linestyle=styles[i],color=cm.Reds(0.9),linewidth=2,label='-box szz')


	#plot differential stress versus axial and radial strain
	print('this is the exx length' + str(len(exx)))
	print('this is the box length' + str(len(box_sdiff)))
	#ax4.plot(-box_avg_exx,-box_avg_diff_stress,linestyle=styles[i],color=cm.gist_rainbow(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box lat strain')
	#ax4.plot(-box_avg_ezz,-box_avg_diff_stress,linestyle=styles[i],color=cm.gist_rainbow(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box ax strain')
	ax4.plot(-box_exx,-box_sdiff,linestyle=styles[i],color=cm.gist_rainbow(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box lat strain')
	ax4.plot(-box_ezz,-box_sdiff,linestyle=styles[i],color=cm.gist_rainbow(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box ax strain')
	#add the volumetric stress vs strain
	#modeled
	ax4.plot(-box_evol,-box_sdiff,linestyle=styles[i],color=cm.Reds(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box vol strain')
	#experimental
	#ax4.plot(-box_ezz,-box_sdiff,linestyle=styles[i],color=cm.gist_rainbow(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box ax strain')
	
 
 
	#plot reaction stress and box average stress versus time 
	ax5.plot(box_time,box_ev,linestyle=styles[i],color=cm.Reds(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' box averaged volumetric stress')
	ax5.plot(box_time,box_p,linestyle=styles[i],color=cm.Blues(colorvec1[colorval1]),linewidth=2,label=testnames[i]+' average box pressure stress')

	ax6.plot(-box_exx,-box_sdiff,linestyle=styles[i],color='black',linewidth=2,label=testnames[i]+' lateral')
	ax6.plot(-box_ezz,-box_sdiff,linestyle=styles[i],color='black',linewidth=2,label=testnames[i]+' axial')

	colorval1 = colorval1+1



################################# PLOTTING EXPERIMENTAL	 DATA ##################################################################
colorvec2=[0.3, 0.5, 0.8]	
colorval2=0



for test in testnames:
	# Open and read the CSV file
	with open(test+'_diffsts_ax_lat_vol.csv', mode='r',encoding='utf-8-sig') as csv_file:  # Replace 'data.csv' with your filename
		reader = csv.reader(csv_file)
		data = list(reader)
		# will need to check individual excel files to see what these indices should be. 
		short_data = data[0:(2900-900)]
	
		# Extract columns as lists
		experimental_stress_ax = [float(row[0])/1000 for row in short_data]  # Column 0
		experimental_strain_ax = [float(row[1]) for row in short_data]  # Column 1
		experimental_strain_lat = [float(row[2]) for row in short_data if len(row) > 2 and row[2]]  # Column 2 (fewer rows)
		experimental_strain_vol = [float(row[3]) for row in short_data if len(row) > 2 and row[2]]  # Column 2 (fewer rows)

		experimental_stress_lat = experimental_stress_ax[:len(experimental_strain_lat)]
		experimental_stress_vol = experimental_stress_ax[:len(experimental_strain_vol)]


		ax4.plot(experimental_strain_ax,experimental_stress_ax,linestyle='--',color=cm.Greys(colorvec2[colorval2]),linewidth=2,label=test+' ax data')
		ax4.plot(experimental_strain_lat,experimental_stress_lat,linestyle='-',color=cm.Greys(colorvec2[colorval2]),linewidth=2,label=test+' lat data')
		ax4.plot(experimental_strain_vol,experimental_stress_vol,linestyle='-',color=cm.Greys(colorvec2[colorval2]),linewidth=2,label=test+' vol data')


################################# Labeling and formatting plots ##################################################################


ax1.set_xlabel(r'time ($\mu$s) ', fontsize=16)
ax1.set_ylabel('F00,F11,F22', fontsize=16)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax1.grid()

ax2.set_xlabel(r'$  pressure (MPa)$', fontsize=16)
ax2.set_ylabel(r' $\sqrt{J_2}$ or   0.57*von Mises stress', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)
ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax2.grid()


ax3.set_xlabel('time', fontsize=16)
ax3.set_ylabel(r'$-\sigma_{xx,yy,zz}$ (GPa)', fontsize=16)
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax3.grid()


ax4.set_xlabel('$-\epsilon_{axial and lateral}$', fontsize=16)
ax4.set_ylabel(r'$-\sigma_{differential}$ (GPa)', fontsize=16)
ax4.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax4.grid()

ax5.set_xlabel('time', fontsize=16)
ax5.set_ylabel(r'$-\sigma_{reaction and box}$ (GPa)', fontsize=16)
ax5.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax5.grid()



ax6.set_xlabel('reaction strain', fontsize=16)
ax6.set_ylabel(r'box stress (GPa)', fontsize=16)
ax6.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax6.grid()


print ('hardening n' + str(strainHardeningn))
# SAVE PLOTS -------------------------------------------------------------------------------


#fig.suptitle('ISR03 Cosine F Table STREN='+str(STREN)+'P='+str(kp)+' I='+str(ki)+ ' D='+str(kd)+' g0='+str(g0)+' g1='+str(g1)+' g2='+str(kd)+' nu='+str(nu),fontsize=20)
#fig.suptitle('ISR03 PEAKI1='+str(PEAKI1)+'P='+str(kp)+' I='+str(ki)+ ' D='+str(kd)+' g0='+str(g0)+' g1='+str(g1)+' g2='+str(kd)+' nu='+str(nu),fontsize=20)
fig.suptitle('ISR03 CR'+str(CR)+ ' k='+str(strainHardeningK)+' n=' +str(strainHardeningn) +'PEAKI1='+str(peakI1)+'STREN='+str(stren)+' \n Fslope='+str(fSlope)+ ' Yslope='+str(ySlope)+' p0='+str(p0)+' b0='+str(b0)+'beta='+str(beta)+'FERR='+str(fractureEnergyReleaseRate)+'frcstss='+str(round(fractureStress,3))+'g0='+str(round(g0,3))+'g1='+str(round(g1,3))+'g2='+str(round(g2,3)),fontsize=20)




fig.tight_layout()
#figname = "ISR03_V1_fracsts:"+str(round(fractureStress,4))+"CR:"+str(CR)+"_beta:"+str(beta)+"fsp:"+str(fSlope)+"ferr:"+str(fractureEnergyReleaseRate)+"hardK:"+str(strainHardeningK)+"hardn:"+str(strainHardeningn)+"stren:"+str(stren)+"g0="+str(round(g0,3))+"g1="+str(round(g1,3))+"g2="+str(round(g2,3))+'new9'
figname='geomechanicsTriaxialCompression_ISR03'
fig.savefig(figname+".png", bbox_inches="tight")
print(figname)
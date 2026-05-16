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
'geoModel_TXC',
]  

labels=[
'geoModel_TXC',
]

styles=[
'-',
'--',
'--',
'-',
'--',
'-',
'--',
'-'
]


# colors=[
# 'blue',
# 'green',
# 'green',
# 'red',
# 'red'
# ]

# input material properties
bulk = 70
shear = 24
C11 = bulk + (4.0/3.0)*shear
EE = 9.0*bulk*shear/(3.0*bulk + shear)
C12 = C11 - 2.0*shear

# initial volume fraction (assumes one material with density rho02)
VF0=np.zeros(len(files))
rho0s=2.7

#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6)) = plt.subplots(3, 2,figsize=(16,16))
fig.suptitle('F table benchmark with uniaxial strain & stress')
#fig2, (ax4,ax5,ax6) = plt.subplots(3, 1,figsize=(16,12))

# plot box sum data and compute initial volume fraction
# for i,file in enumerate(files):
# 	print(file)
# 	reactionFile=runLocation+file+'/boxSumHistory.csv'

# 	# time,sxx,syy,szz,sxy,syz,sxz,rho,energy
# 	data = np.genfromtxt(reactionFile, delimiter=',')
# 	time = data[:,0]
# 	sxx = data[:,1]
# 	syy = data[:,2]
# 	szz = data[:,3]
# 	sxy = data[:,4]
# 	syz = data[:,5]
# 	sxz = data[:,6]
# 	rho = data[:,7]
# 	e = data[:,8]

# 	# this hids all non-monotic time entries, so restart data files are cleaned:
# 	maxt = 0.0
# 	mask = np.ones(len(time), dtype=bool)
# 	for ii,t in enumerate(time):
# 		if (t<=maxt):
# 			mask[ii] = False
# 		else:
# 			maxt = t
# 	time = time[mask,...]
# 	sxx = sxx[mask,...]
# 	syy = syy[mask,...]
# 	szz = szz[mask,...]
# 	sxy = sxy[mask,...]
# 	syz = syz[mask,...]
# 	sxz = sxz[mask,...]
# 	rho = rho[mask,...]
# 	e = e[mask,...]


# 	# p=(-1.0/3.0)*(sxx+syy+szz)
# 	# vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )
# 	ax3.plot(time,1000*sxx,linestyle=':',color=cm.gist_rainbow(0),linewidth=1,label='box: sxx')
# 	ax3.plot(time,1000*syy,linestyle=':',color=cm.gist_rainbow(0.33),linewidth=1,label='box: syy')
# 	ax3.plot(time,1000*szz,linestyle=':',color=cm.gist_rainbow(0.67),linewidth=1,label='box: szz')


# print("VF0 = ",VF0) 


for i,file in enumerate(files):
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

	p=(-1.0/3.0)*(sxx+syy+szz)
	rho=rho0/F00

	vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )


    # stress is scaled from simulation units (GPa) to plotting units (MPa)
	ax1.plot(time,exx,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='exx')
	ax1.plot(time,eyy,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=1,label='eyy')
	ax1.plot(time,ezz,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='ezz')


	#ax2.plot(p,vm,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='reaction')

	# ax2.plot(exx,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	# ax2.plot(exx,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	# ax2.plot(exx,-1000*Rym/Ay,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=2,label='y_minus')
	# ax2.plot(exx,1000.0*Ryp/Ay,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=2,label='y_plus')
	# ax2.plot(exx,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='z_minus')
	# ax2.plot(exx,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')

	# ax2.plot(exx,1e3*C11*exx,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Constrained (C11) Modulus')
	# ax2.plot(exx,1e3*EE*exx,linestyle='--',color=cm.gnuplot2(0),linewidth=1,label='Unconstrained (Young\'s) Modulus')
	# ax2.plot(exx,1e3*C12*exx,linestyle='-.',color=cm.gnuplot2(0),linewidth=1,label='C12')

	# ax3.plot(time,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	# ax3.plot(time,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	# ax3.plot(time,-1000*Rym/Ay,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=2,label='y_minus')
	# ax3.plot(time,1000.0*Ryp/Ay,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=2,label='y_plus')
	# ax3.plot(time,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=2,label='z_minus')
	# ax3.plot(time,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')
 

    # the strain are plotted relative to the hydrostatic state at t=tStop/2 (change this if stress table time changes)
	i0 = int(len(ezz)/2)
 
	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)
	stren = job.STREN
	peakI1 = job.PEAKI1
	ySlope = job.YSLOPE
	fSlope = job.FSLOPE
	fSlopeFailed = job.FSLOPEFAILED
	strainHardeningK = job.strainHardeningK
	fractureEnergyReleaseRate = job.fractureEnergyReleaseRate
	fractureStress = job.fractureStress
 
 
	print('[stren, peakI1, ySlope,fSlope,fSlopeFailed, Gf, sigmaF] = ',stren, peakI1, ySlope,fSlope,fSlopeFailed,fractureEnergyReleaseRate,fractureStress)
 
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

 
	print('a1,a2,a3,a4 = ',a1,a2,a3,a4)
 
	I1plot = np.linspace(p0,peakI1,100)
	pPlot = np.array([ -I1/3. for I1 in I1plot ])

	Ff_MPa = 1000.*np.array([Ff(I1,a1,a2,a3,a4) for I1 in I1plot ])
	FfFc_MPa = Ff_MPa*np.array([Fc(I1,p0,peakI1,CR) for I1 in I1plot ])
	ax4.plot( 1000.*pPlot,FfFc_MPa,linestyle='--',color='b',linewidth=3,label='Fc*Ff')
 
	I1Fplot = np.linspace(p0,0,100)
	pFPlot = np.array([ -I1/3. for I1 in I1Fplot ])
	FfFailed_MPa = 1000.*np.array([Ff(I1,a1f,a2f,a3f,a4f) for I1 in I1Fplot ])
	FfFailedFc_MPa = FfFailed_MPa*np.array([Fc(I1,p0,0,CR) for I1 in I1Fplot ])
	ax4.plot( 1000.*pFPlot,FfFailedFc_MPa,linestyle='--',color='r',linewidth=3,label='Fc*Ff_failed')	

	I1Hplot = np.linspace(p0,peakI1h,100)
	pHPlot = np.array([ -I1/3. for I1 in I1Hplot ])
	FfHardened_MPa = 1000.*np.array([Ff(I1,a1h,a2h,a3h,a4h) for I1 in I1Hplot ])
	FfHardenedFc_MPa = FfHardened_MPa*np.array([Fc(I1,p0,peakI1h,CR) for I1 in I1Hplot ])	
	ax4.plot( 1000.*pHPlot,FfHardenedFc_MPa,linestyle='--',color='c',linewidth=3,label='Fc*Ff_hardened')	
 
  
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
    
		box_p = (-1.0/3.0)*(box_sxx+box_syy+box_szz)
		box_vm = np.sqrt(0.5*( np.square(box_sxx-box_syy)+np.square(box_syy-box_szz)+np.square(box_szz-box_sxx) ) )
		box_ev = -np.log(box_rho[0] / np.array(box_rho))
		box_evp = box_epxx + box_epyy + box_epzz

	#print('box data',box_time, box_sxx)  
        # plot strain history
		ax1.plot(box_time,box_evp,linestyle='-',color=cm.gist_rainbow(.6),linewidth=3,label='box evp')

        # plot plastic strain and damage time history
		ax2.plot(box_time,box_epxx,linestyle='-',color=cm.gist_rainbow(0),linewidth=3,label='box ep_xx')		
		ax2.plot(box_time,box_epyy,linestyle='-',color=cm.gist_rainbow(.1),linewidth=3,label='box ep_yy')
		ax2.plot(box_time,box_epzz,linestyle='-',color=cm.gist_rainbow(.2),linewidth=3,label='box ep_zz')
		# ax2.plot(box_time,box_epxy,linestyle='-',color=cm.gist_rainbow(.3),linewidth=3,label='box ep_xy')
		# ax2.plot(box_time,box_epyz,linestyle='-',color=cm.gist_rainbow(.4),linewidth=3,label='box ep_yz')
		# ax2.plot(box_time,box_epxz,linestyle='-',color=cm.gist_rainbow(.5),linewidth=3,label='box ep_xz')
		ax2.plot(box_time,box_evp,linestyle='-',color=cm.gist_rainbow(.6),linewidth=3,label='box evp')		
		
		
        # plot stress history, make sure lateral stress is keeping constant
		ax3.plot(box_time,1000.*box_sxx,linestyle='-',color=cm.gist_rainbow(0),linewidth=3,label='box sxx')		
		ax3.plot(box_time,1000.*box_syy,linestyle='-',color=cm.gist_rainbow(.2),linewidth=3,label='box syy')		
		ax3.plot(box_time,1000.*box_szz,linestyle='-',color=cm.gist_rainbow(.4),linewidth=3,label='box szz')		
		ax3.plot(box_time,1000.*box_p,linestyle='-',color=cm.gist_rainbow(.6),linewidth=3,label='box p')		
		ax3.plot(box_time,1000.*box_vm,linestyle='-',color=cm.gist_rainbow(.8),linewidth=3,label='box vm')		
  
		box_I1 = np.array([ -3.*p for p in box_p ])
		box_Ff_MPa = 1000.*np.array([ a1 - a3*np.exp(a2*I1) - a4*I1 for I1 in box_I1 ])
		ax3.plot(box_time,np.sqrt(3.)*box_Ff_MPa,linestyle='--',color='k',linewidth=2,label='box Ff')
    
        # plt pressure-shear, it would be good to add the yield surface to this plot.
		ax4.plot(1000.*box_p,np.sqrt(1/3)*1000.*box_vm,linestyle='-',color=cm.gist_rainbow(5),linewidth=3,label='box')
		ax5.plot(box_time,box_damage,linestyle='-',color='r',linewidth=3,label='box damage')		

		ax6.plot(-100.*box_ezz,-1000.*(box_szz-box_sxx),linestyle='-',color=cm.gist_rainbow(0.27),linewidth=1,label='-ezz')
		ax6.plot(-100.*box_exx,-1000.*(box_szz-box_sxx),linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='-exx')
		ax6.plot(-100.*(box_exx+box_eyy+box_ezz),-1000.*(box_szz-box_sxx),linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='-ev')
  
		print('Percent complete: ',box_time[-1]/job.stopTime)

ax1.set_xlabel('time')
#ax1.set_xlim(0.0,-1.05*min(exx[1:]))
ax1.set_ylabel('Total Strain')
#ax1.set_yscale("log")
#ax1.set_ylim(-2.5,2.5)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax1.grid()

ax2.set_xlabel('time')
ax2.set_xlabel('Plastic strain')
ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax2.grid()

ax3.set_xlabel('time')
ax3.set_xlabel('Stress (MPa)')
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax3.grid()

ax4.set_xlabel('Pressure (MPa)')
ax4.set_ylabel('rootJ2 shear stress = sqrt(1/3)*von Mises (MPa)')
ax4.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax4.grid()

ax5.set_xlabel('time', fontsize=16)
ax5.set_ylabel('damage', fontsize=16)
ax5.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')

ax6.set_xlabel('lateral strain | axial strain (%)', fontsize=16)
ax6.set_ylabel('differential stress (MPa)', fontsize=16)
ax6.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')


fig.tight_layout()
#plt.show()
fig.savefig("geoModel_TXC.png", bbox_inches="tight")
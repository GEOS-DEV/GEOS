# -*- coding: utf-8 -*-
"""
Created on Wed Mar 1 09:00:00 2017
@author: homel1
Geometry object functions for the particle file writer.
"""
import numpy as np                   # math stuff
import matplotlib.pyplot as plt

import matplotlib.cm as cm

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
# 'geoModel_noBuckling_unconfinedCompression2D',
# 'geoModel_isotropicBuckling_unconfinedCompression2D',
# 'geoModel_anisotropicBuckling_unconfinedCompression2D',
'geoModel_noBuckling_unconfinedCompression2DQS',
'geoModel_isotropicBuckling_unconfinedCompression2DQS',
'geoModel_anisotropicBuckling_unconfinedCompression2DQS'
]  

# labels=[
# r'no buckling: $\dot{\varepsilon}=5.\times10^4/s$',
# r'isotropic buckling: $\dot{\varepsilon}=5.\times10^4/s$',
# r'anisotropic buckling [101]: $\dot{\varepsilon}=5.\times10^4/s$',
# r'no buckling: $\dot{\varepsilon}=5.\times10^3/s$',
# r'isotropic buckling: $\dot{\varepsilon}=5.\times10^3/s$',
# r'anisotropic buckling [101]: $\dot{\varepsilon}=5.\times10^3/s$'
# ]

labels=[
# r'none: $\dot{\varepsilon}=5.\times10^4/s$',
# r'iso: $\dot{\varepsilon}=5.\times10^4/s$',
# r'aniso [101]: $\dot{\varepsilon}=5.\times10^4/s$',
r'none: $\dot{\varepsilon}=5.\times10^3/s$',
r'iso: $\dot{\varepsilon}=5.\times10^3/s$',
r'aniso [101]: $\dot{\varepsilon}=5.\times10^3/s$'
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
fig, axes = plt.subplots(2, 2,figsize=(15,12))
fig.suptitle('Unconfined Plane-Strain Compression')
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
	color = cm.gist_rainbow(i/len(files))
	label = labels[i]
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
	# axes[0][0].plot(time,F00,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='F00')
	axes[0][0].plot(time,F11,linestyle='-',color=color,linewidth=1,label=r'$F_{11}$')
	#axes[0][0].plot(time,F22,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='F22')

	#axes[1][0].plot(p,vm,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='reaction')

	# ax2.plot(exx,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	# ax2.plot(exx,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	# ax2.plot(exx,-1000*Rym/Ay,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=2,label='y_minus')
	# ax2.plot(exx,1000.0*Ryp/Ay,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=2,label='y_plus')
	# ax2.plot(exx,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='z_minus')
	# ax2.plot(exx,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')

	# ax2.plot(exx,1e3*C11*exx,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Constrained (C11) Modulus')
	# ax2.plot(exx,1e3*EE*exx,linestyle='--',color=cm.gnuplot2(0),linewidth=1,label='Unconstrained (Young\'s) Modulus')
	# ax2.plot(exx,1e3*C12*exx,linestyle='-.',color=cm.gnuplot2(0),linewidth=1,label='C12')

	# axes[0][1].plot(time,-1000*Rxm/Ax,linestyle='-',color=cm.gist_rainbow(0),linewidth=2,label='x_minus')
	# axes[0][1].plot(time,1000.0*Rxp/Ax,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=2,label='x_plus')
	axes[0][1].plot(time,1000*Rym/Ay,linestyle='-',color=color,linewidth=2,label=r'$\sigma_y^-$')
	axes[0][1].plot(time,-1000.0*Ryp/Ay,linestyle='--',color=lighten_color(color,1.5),linewidth=2,label=r'$\sigma_y^+$')
	#axes[0][1].plot(time,-1000*Rzm/Az,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=2,label='z_minus')
	#axes[0][1].plot(time,1000.0*Rzp/Az,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=2,label='z_plus')
 
 #print(file)
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
		box_epyx = data[1:,14]
		box_epxz = data[1:,15]
		box_epxy = data[1:,16]
		box_volFrac = data[1:,17]
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
		box_sxy = box_sxy[mask,...]
		box_syz = box_syz[mask,...]
		box_sxz = box_sxz[mask,...]
		box_rho = box_rho[mask,...]
		box_e = box_e[mask,...]
		box_damage = box_damage[mask,...]
		box_volFrac = box_volFrac[mask,...]
		box_p = (-1.0/3.0)*(box_sxx+box_syy+box_szz)
		box_vm=np.sqrt(0.5*( np.square(box_sxx-box_syy)+np.square(box_syy-box_szz)+np.square(box_szz-box_sxx) ) )

	#print('box data',box_time, box_sxx)  
		#axes[0][1].plot(box_time,1000.*box_sxx,linestyle='-',color=cm.gist_rainbow(0),linewidth=3,label='box sxx')		
		axes[0][1].plot(box_time,-1000.*box_syy,linestyle='-',color=color,linewidth=3,label=r'$\overline{\sigma_y}$ '+labels[i])		
		#axes[0][1].plot(box_time,1000.*box_szz,linestyle='-',color=cm.gist_rainbow(.4),linewidth=3,label='box szz')		
		
		#axes[0][1].plot(box_time,1000.*box_p,linestyle='-',color=cm.gist_rainbow(.6),linewidth=3,label='box p')		
		#axes[0][1].plot(box_time,1000.*box_vm,linestyle='-',color=cm.gist_rainbow(.8),linewidth=3,label='box vm')		
		axes[1][0].plot(-eyy,0.8*1000*Rym/Ay,linestyle='-',color=color,linewidth=2,label=r'$\sigma_y^-$')
		axes[1][0].plot(-eyy,-0.8*1000*Ryp/Ay,linestyle='--',color=lighten_color(color,1.5),linewidth=2,label=r'$\sigma_y^+$')
  
		axes[1][1].plot(box_rho,box_p,linestyle='-',color=color,linewidth=3,label=r'$p$')
		axes[1][1].plot(box_rho,box_vm,linestyle='--',color=color,linewidth=3,label=r'$\tau$')



# von Mises vs density
axes[0][0].set_xlabel('time')
#ax1.set_xlim(0.0,-1.05*min(exx[1:]))
axes[0][0].set_ylabel('F')
#ax1.set_yscale("log")
#ax1.set_ylim(0.001,5.0)
#ax1.set_ylim(-2.5,2.5)
axes[0][0].legend( loc="best",fontsize='x-small')
axes[0][0].grid()

# plot solid eos
# js = np.arange(0.8, 1.0001, 0.0001)
# ps = -bulk*np.log(js)
# ax2.plot(1/js, ps, 'k-',linewidth=2,label='Solid EOS')

axes[1][0].set_xlabel('Compressive Strain', fontsize=16)
#ax2.set_xlim(0,1.2)
#ax2.set_ylim(0,600)
#ax2.set_yscale("log")
#ax2.set_ylim(0.001,1000)

axes[1][0].set_ylabel(r'$\sigma_{vm}$ (MPa)', fontsize=16)
axes[1][0].tick_params(axis='both', which='major', labelsize=16)
axes[1][0].tick_params(axis='both', which='minor', labelsize=16)

axes[1][0].legend(loc="upper right",fontsize='medium')
axes[1][0].grid()


axes[0][1].set_xlabel('time', fontsize=16)
axes[0][1].set_ylabel(r'$\sigma$ (MPa)', fontsize=16)
axes[0][1].legend(loc="best",fontsize='medium')


axes[1][1].set_xlabel('density', fontsize=16)
axes[1][1].set_ylabel(r'$\sigma$ (MPa)', fontsize=16)
axes[1][1].legend(loc="best",fontsize='medium')



# ax3.set_xlim(0.00,25)
# ax3.set_ylim(0.00,200)
#ax3.set_yscale("log")
#ax3.set_ylim(0.001,1000)


fig.tight_layout()
#plt.show()
fig.savefig("geoModel_UnconfinedBuckling.png", bbox_inches="tight")
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


# parent directory of file path:
runLocation='/p/lustre1/crook5/geosxRuns/'

files=[
'ceramicTXC'
]  

labels=[
'ceramic'
]

area = [
1.0
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


# colors=[
# 'blue',
# 'green',
# 'green',
# 'red',
# 'red'
# ]

#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, (ax1, ax2, ax3) = plt.subplots(3, 1,figsize=(12,9))
fig.suptitle('Uniaxial strain verification of von Mises Model')
#fig2, (ax4,ax5,ax6) = plt.subplots(3, 1,figsize=(16,12))

# plot box sum data and compute initial volume fraction
for i,file in enumerate(files):
	print(file)
	reactionFile=runLocation+file+'/boxAverageHistory.csv'

	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)
	bulk = job.bulk
	shear = job.shear
	rho0 = job.density
	compressiveStrength = job.compressiveStrength
	tensileStrength = job.tensileStrength
	maximumStrength = job.maximumStrength
	p0 = job.p0

	print('bulk = ',bulk,", shear = ",shear,", rho0 = ",rho0,", compressiveStrength = ",compressiveStrength)

	# time,sxx,syy,szz,sxy,syz,sxz,rho,energy
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[:,0]
	sxx = data[:,1]
	syy = data[:,2]
	szz = data[:,3]
	sxy = data[:,4]
	syz = data[:,5]
	sxz = data[:,6]
	rho = data[:,7]
	e = data[:,8]

	# this hids all non-monotic time entries, so restart data files are cleaned:
	maxt = 0.0
	mask = np.ones(len(time), dtype=bool)
	for ii,t in enumerate(time):
		if (t<=maxt or np.isnan(t)):
			mask[ii] = False
		else:
			maxt = t
	time = time[mask,...]
	sxx = sxx[mask,...]
	syy = syy[mask,...]
	szz = szz[mask,...]
	sxy = sxy[mask,...]
	syz = syz[mask,...]
	sxz = sxz[mask,...]
	rho = rho[mask,...]
	e = e[mask,...]

	J = rho0*np.reciprocal(rho)

	#ev = np.log(J)
	mu = 0.5773502691896258
	Yt = job.tensileStrength
	Yc = job.compressiveStrength
	Ymax = job.maximumStrength

	Yt0 = max( 0.5*Yt, min( 2.0*Yt, (3.*Yc*Yt)/(2.*Yc + Yt) ) )
	print(Yt0, 0.5*Yt, 2.0*Yt, (3.*Yc*Yt)/(2.*Yc + Yt))
	Yt0 = min( Yt0, ( 3.0 * Yc - Yc * mu ) / ( 3 + mu ) )

	pmin0 = (-2*Yc*Yt0)/(3.*(Yc - Yt0))

	x1 = Yc/3
	x2 = Ymax/mu
	y1 = Yc
	y2 = Ymax
	m1 = (Yc - Yt0)/(Yc/3+Yt0/3)
	beta = (m1*x1 - m1*x2 - y1 + y2)/(y1 - y2)

	# max pressure in simulation.

	p=(-1.0/3.0)*(sxx+syy+szz)
	ax1.plot(-np.log(J),-bulk*np.log(J),linestyle='-',color='y',linewidth=4,label='-k ln(J)')
	ax1.plot(-np.log(J),p0*np.ones(J.size),linestyle='--',color='r',linewidth=2,label='p0')
	ax1.plot(-np.log(J),p,linestyle='--',color='k',linewidth=2,label='Sim output: p')

	pmax = np.max(p)


	# plot initial yield surface in TXC: Y(P), D=0
	prange = np.linspace(pmin0, pmax, num=100)
	Y=maximumStrength*np.ones(prange.size)

		# plot ymax
	ax2.plot(prange,Y,linestyle='--',color='r',linewidth=2,label='maximumStrength')

		## plot Y(p) d=1
	for ip, pp in enumerate(prange):
		if (pp<Yc/3):
			Y[ip] = m1*(pp-pmin0)
		elif (pp < Ymax/mu):
			Y[ip] = (m1*(pp - x2)*((pp - x2)/(x1 - x2))**beta + m1*(-x1 + x2) + (1 + beta)*y1)/(1 + beta)
		else:
			Y[ip] = Ymax

	ax2.plot(prange,Y,linestyle='-',color='y',linewidth=2,label=r'$Y(p)$, $D=0$')




	prange = np.linspace(0, pmax, num=100)
	Y=maximumStrength*np.ones(prange.size)
	## plot Y(p) d=1
	for ip, pp in enumerate(prange):
		Y[ip] = min(mu*pp,maximumStrength)

	ax2.plot(prange,Y,linestyle='-',color='b',linewidth=2,label=r'$Y(p)$, $D=1$')


	#ax2.plot(psol,vmsol,linestyle='-',color='g',linewidth=4,label='Solution')
	vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )
	ax2.plot(p,vm,linestyle='--',color='k',linewidth=2,label='Sim output: vm stress')

	pd = np.array([-Yt0/3., -Yt/3., Yc/3.])
	yd = np.array([Yt0, Yt, Yc])
	ax2.plot(pd, yd, 'o', color='black');

	# ax1.plot(-np.log(J),-bulk*np.log(J),linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='-k ln(J)')
	
	ax3.plot(time,sxx,linestyle=':',color=cm.gist_rainbow(0),linewidth=2,label='Sim output: sxx')
	ax3.plot(time,syy,linestyle=':',color=cm.gist_rainbow(0.2),linewidth=2,label='Sim output: syy')
	ax3.plot(time,szz,linestyle=':',color=cm.gist_rainbow(0.4),linewidth=2,label='Sim output: szz')
	ax3.plot(time,vm,linestyle=':',color=cm.gist_rainbow(0.6),linewidth=2,label='Sim output: mises')
	ax3.plot(time,-p,linestyle=':',color=cm.gist_rainbow(0.8),linewidth=2,label='Sim output: mean')

# print("VF0 = ",VF0) 


for i,file in enumerate(files):

	print(file)

	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)
	bulk = job.bulk
	shear = job.shear
	rho0 = job.density
	compressiveStrength = job.compressiveStrength

	reactionFile=runLocation+file+'/reactionHistory.csv'

	# time,F00,F11,F22,lx,ly,lz,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
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
	maxt = 0.0
	mask = np.ones(len(time), dtype=bool)
	for ii,t in enumerate(time):
		if (t<=maxt):
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

	#vm=np.sqrt(0.5*( np.square(sxx-syy)+np.square(syy-szz)+np.square(szz-sxx) ) )


    # stress is scaled from simulation units (GPa) to plotting units (MPa)
	# ax1.plot(time,F00,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='F00')
	# ax1.plot(time,F11,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=1,label='F11')
	# ax1.plot(time,F22,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='F22')

	#ax2.plot(time,vm,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='VM (reactions)')

	# ax2.plot(exx,-1000*Rxm/A0,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='x_minus')
	# ax2.plot(exx,1000.0*Rxp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=1,label='x_plus')
	# ax2.plot(exx,-1000*Rym/A0,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=1,label='y_minus')
	# ax2.plot(exx,1000.0*Ryp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=1,label='y_plus')
	# ax2.plot(exx,-1000*Rzm/A0,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='z_minus')
	# ax2.plot(exx,1000.0*Rzp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=1,label='z_plus')

	# # ax2.plot(exx,1e3*C11*exx,linestyle='-',color=cm.gnuplot2(0),linewidth=1,label='Constrained (C11) Modulus')
	# # ax2.plot(exx,1e3*EE*exx,linestyle='--',color=cm.gnuplot2(0),linewidth=1,label='Unconstrained (Young\'s) Modulus')
	# # ax2.plot(exx,1e3*C12*exx,linestyle='-.',color=cm.gnuplot2(0),linewidth=1,label='C12')

	# ax3.plot(time,-1000*Rxm/A0,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='x_minus')
	# ax3.plot(time,1000.0*Rxp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0),1.5),linewidth=1,label='x_plus')
	# ax3.plot(time,-1000*Rym/A0,linestyle='-',color=cm.gist_rainbow(0.33),linewidth=1,label='y_minus')
	# ax3.plot(time,1000.0*Ryp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0.33),1.5),linewidth=1,label='y_plus')
	# ax3.plot(time,-1000*Rzm/A0,linestyle='-',color=cm.gist_rainbow(0.67),linewidth=1,label='z_minus')
	# ax3.plot(time,1000.0*Rzp/A0,linestyle='--',color=lighten_color(cm.gist_rainbow(0.67),1.5),linewidth=1,label='z_plus')



# von Mises vs density
ax1.set_xlabel('-ln*(J)', fontsize=16)
#ax1.set_xlim(0.0,-1.05*min(exx[1:]))
ax1.set_ylabel('p, (GPa)', fontsize=16)
#ax1.set_yscale("log")
#ax1.set_ylim(0.001,5.0)
#ax1.set_ylim(-2.5,2.5)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax1.grid()

# plot solid eos
# js = np.arange(0.8, 1.0001, 0.0001)
# ps = -bulk*np.log(js)
# ax2.plot(1/js, ps, 'k-',linewidth=2,label='Solid EOS')

ax2.set_xlabel('pressure (GPa)', fontsize=16)
#ax2.set_xlim(0,1.2)
#ax2.set_ylim(0,600)
#ax2.set_yscale("log")
#ax2.set_ylim(0.001,1000)

ax2.set_ylabel(r'$\sigma_{vm}$ (GPa)', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)

ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax2.grid()


ax3.set_xlabel('time', fontsize=16)
ax3.set_ylabel(r'$\sigma$ (GPa)', fontsize=16)
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')

#ax3.set_yscale("log")
#ax3.set_ylim(0.001,1000)


fig.tight_layout()
#plt.show()
fig.savefig("ceramicTXC.png", bbox_inches="tight")
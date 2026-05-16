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
import csv



#please enter the name of the exerimental test data you're trying to constrain (TP1, TP2, TP3)
experiments = ['TP1','TP2','TP3']
lat_diff_stress_experiments = [39.39,25.02,39.32]
lat_diff_stress_experiments = [item/1000 for item in lat_diff_stress_experiments]
ax_diff_stress_experiments = [39.34,24.85,39.27]
ax_diff_stress_experiments = [item/1000 for item in ax_diff_stress_experiments]

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
'geomechanicsTriaxialCompression_TP1',
'geomechanicsTriaxialCompression_TP2',
'geomechanicsTriaxialCompression_TP3'
]  

labels=[
'TP1',
'TP2',
'TP3'
]
# Material model parameters.  Will be read from input file for plotting.
#density = 2.648
bulk = 36.3
shear = 26.0
sigmaYield = 2.27
sigmaFail = 0.449
sigmaMax = 5.0
crackSpeed = 1.8

youngsModulus = (9.*bulk*shear)/(3.*bulk+shear)
strainToFailure = sigmaFail/youngsModulus
stretch = np.exp(2.0*strainToFailure) # stretch to 50% past the strain to failure.

# useDamageAsInitialStrength must be enabled for the painted damage field to be used as strengthScale 
useDamageAsInitialStrength=0               
weibullSeed=1
weibullModulus=6.0
weibullVolume=1.0
flawSize=0.2

sampleHeight = 1. # mm
sampleWidth = 1.  # mm



styles=[
'-',
'--',
'-.',
'-',
'--',
'-.',
'-',
'--',
'-.'
]



#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]
fig, ((ax1, ax2), (ax3, ax4), (ax5,ax6)) = plt.subplots(3, 2,figsize=(24,24))
fig.suptitle('triaxial compaction, strain=0.05 stress=20.7 kd=0.01, stoptime = 10000',fontsize='large')

# BOX AVERAGE HISTORY FILE -------------------------------------------------------------------------------
# plot box sum data and compute initial volume fraction
for i,file in enumerate(files):
	print(file)
	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)

	density = job.density 
	b0 = job.b0 
	b1 = job.b1 
	b2 = job.b2 
	b3 = job.b3 
	b4 = job.b4 
	g0 = job.g0 
	g1 = job.g1 
	g2 = job.g2 
	g3 = job.g3 
	g4 = job.g4 
	PEAKI1 = job.PEAKI1 
	FSLOPE = job.FSLOPE 
	STREN = job.STREN 
	YSLOPE = job.YSLOPE 
	BETA_nonassociativity = job.BETA_nonassociativity 
	p0 = job.p0 
	p1 = job.p1 
	p2 = job.p2 
	p3 = job.p3 
	CR = job.CR 
	Kf = job.Kf 
	pfi = job.pfi 
	strainHardeningK = job.strainHardeningK
	strainHardeningn = job.strainHardeningn
	A = job.creepA
	B = job.creepB
	C = job.creepC
	creepc0 = job.creepc0
	creepc1 = job.creepc1
	confiningPressure = job.confiningPressure
	stopTime = job.stopTime

	reactionFile=runLocation+file+'/boxAverageHistory.csv'


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
 
 

	#print(type(box_avg_Sxx))
 
	
	#print('box diff stress' +str(box_avg_diff_stress))

	#print('this is box_avg_exx'+str(box_avg_exx))
	#print('this is box_avg_vol strain'+str(box_avg_vol_strain))

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
	reactionFile=runLocation+file+'/reactionHistory.csv'
	# area (mm^2)
	A0=sampleWidth**2

	# time,F00,F11,F22,Rx-,Rx+,Ry-,Ry+,Rz-,Rz+
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[1:,0]
	F00 = data[1:,1]
	F11 = data[1:,2]
	F22 = data[1:,3]
	#print(str(F11))
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
	#ev = -np.log(F00*F11*F22)
	ev = (-1. )*(exx + eyy + ezz)
	#this is how they calculate it in kikibas et al 2022. So this it the volumetric strain in the plot
 
	# The contact areas
	#Ax=A0*F11*F22
	#Ay=A0*F00*F22
	#Az=A0*F00*F11

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


	#volumetric strain
	#ev = (1.0)*(exx+eyy+ezz)

	#vm = np.sqrt( 0.5*( np.square(sxx-syy) + np.square(syy-szz) + np.square(szz-sxx) ) )


    # stress is scaled from simulation units (GPa) to plotting units (MPa)
	ax1.plot(time,F00,linestyle='-',color=cm.gist_rainbow(0),linewidth=1,label='F00')
	ax1.plot(time,F11,linestyle='-',color=cm.gist_rainbow(0.25),linewidth=1,label='F11')
	ax1.plot(time,F22,linestyle='-',color=cm.gist_rainbow(0.75),linewidth=1,label='F22')


# PLOTTING SIMULATED DATA -------------------------------------------------------------------------------


	# Plot Yield Surface (non-hardening)
	a1 = STREN;
	a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
	a3 = (STREN-YSLOPE*PEAKI1)*np.exp(-a2*PEAKI1);
	a4 = YSLOPE;

	X = p0
	I1Plot = np.linspace(X, PEAKI1, num=100)
	Ff = np.sqrt(3.)*(a1 - a3*np.exp(a2*I1Plot) - a4*I1Plot)
	Kappa  = PEAKI1-CR*(PEAKI1-X)
	Fc = np.array([ (1 if I1 > Kappa else ( 0 if I1 < X else np.sqrt(1.0 - ((Kappa-I1)/(Kappa-X))*((Kappa-I1)/(Kappa-X))) ) ) for I1 in I1Plot])
	ax2.plot((-1./3.)*I1Plot,Ff,'--k',linewidth=2,label='Ff (non-hardening)')
	ax2.plot((-1./3.)*I1Plot,Fc*Ff,'--b',linewidth=2,label='Ff*Fc (non-hardening)')



	# Plot Limit Surface (hardening)
	a1 = STREN + strainHardeningK;
	a2 = (FSLOPE-YSLOPE)/(STREN-YSLOPE*PEAKI1);
	a3 = (STREN-YSLOPE*PEAKI1)*np.exp(-a2*PEAKI1);
	a4 = YSLOPE;
	Ff = np.sqrt(3.)*(a1 - a3*np.exp(a2*I1Plot) - a4*I1Plot)
 
	#ax2.plot((-1./3.)*I1Plot,Ff,'-k',linewidth=2,label='Ff (hardening)')

	# plot model results
	ax2.plot(box_avg_vol_stress,box_avg_diff_stress,linestyle=styles[i],color=colors[i],linewidth=2,label='model')
 
	#plot data from experiments
	ax2.plot(confiningPressure,((lat_diff_stress_experiments[i])/np.sqrt(3)),'o',color=colors[i],label='lateral experimental value')
	ax2.plot(confiningPressure,((ax_diff_stress_experiments[i])/np.sqrt(3)),'s',color=colors[i],label='axial experimental value')

	ax3.plot(box_time,-box_avg_Szz,linestyle=styles[i],color=cm.gist_rainbow(0.2),linewidth=2,label='-Bsxx')
	ax3.plot(box_time,-box_avg_Syy,linestyle=styles[i],color=cm.gist_rainbow(0.4),linewidth=2,label='-Bsyy')
	ax3.plot(box_time,-box_avg_Szz,linestyle=styles[i],color=cm.gist_rainbow(0.6),linewidth=2,label='-Bszz')

	#plot differential stress versus axial and radial strain
	ax4.plot(-box_exx,box_avg_diff_stress,linestyle=':',color=colors[i],linewidth=2,label=labels[i]+' model lat strain')
	ax4.plot(-box_ezz,box_avg_diff_stress,linestyle='--',color=colors[i],linewidth=2,label=labels[i]+'model ax strain')
	#ax4.plot(time,-szz,linestyle=styles[i],color=cm.gist_rainbow(0.6),linewidth=2,label='-szz')
 	
	#plot reaction stress and box average stress versus time 
	ax5.plot(box_time,box_avg_vol_stress,linestyle=styles[i],color=cm.gist_rainbow(0.2),linewidth=2,label='box averaged volumetric stress')
	
	#ax6.plot(time,-sxx,linestyle=styles[i],color=cm.gist_rainbow(0.2),linewidth=2,label='-sxx')
	#ax6.plot(time,-syy,linestyle=styles[i],color=cm.gist_rainbow(0.4),linewidth=2,label='-syy')
	#ax6.plot(time,-szz,linestyle=styles[i],color=cm.gist_rainbow(0.6),linewidth=2,label='-szz')

# PLOTTING EXPERIMENTAL	 DATA -------------------------------------------------------------------------------
rainbowvec=[0.0, 0.2, 0.4]	
rainbowval=0

for i,test in enumerate(experiments):
		experimental_stress = []
		experimental_ax_strain = []
		with open(test+'_diff_stress_vs_ax_strain_2024.csv', mode='r') as csv_file:
			csv_reader=csv.reader(csv_file)
			next(csv_reader,None)
			for row in csv_reader:
				experimental_ax_strain.append(float(row[0]))
				experimental_stress.append(float(row[1])/1000)
		ax4.plot(experimental_ax_strain,experimental_stress,linestyle='-',color=colors[i],linewidth=2,label=test+' ax data')

		experimental_stress = []
		experimental_lat_strain = []
		with open(test+'_diff_stress_vs_lat_strain_2024.csv', mode='r') as csv_file:
			csv_reader=csv.reader(csv_file)
			next(csv_reader,None)
			for row in csv_reader:
				experimental_lat_strain.append(float(row[0]))
				experimental_stress.append(float(row[1])/1000)
		ax4.plot(experimental_lat_strain,experimental_stress,linestyle='-',color=colors[i],linewidth=2,label=test+' lat data')
		rainbowval=rainbowval+1


# LABELING/FORMATTING PLOTS -------------------------------------------------------------------------------

ax1.set_xlabel(r'time ($\mu$s) ', fontsize=16)
#ax1.set_xlim(0.0,-1.05*min(exx[1:]))
ax1.set_ylabel('F00,F11,F22', fontsize=16)
#ax1.set_yscale("log")
#ax1.set_ylim(0.0001,1.0)
#ax1.set_ylim(-2.5,2.5)
ax1.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='x-small')
ax1.grid()

ax2.set_xlabel(r'$pressure$', fontsize=16)
#ax2.set_xlim(0,1.2)
#ax2.set_ylim(0.,1.5*sigmaYield)
#ax2.set_yscale("log")
#ax2.set_ylim(0.001,1000)

ax2.set_ylabel(r'$von Mises stress$ (GPa)', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)
ax2.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
#ax2.grid()


ax3.set_xlabel('time', fontsize=16)
ax3.set_ylabel(r'$-\sigma_{xx,yy,zz}$ (GPa)', fontsize=16)
ax3.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax3.grid()
#ax3.set_yscale("log")
#ax3.set_ylim(0.,0.05)

ax4.set_xlabel('$-\epsilon_{axial and lateral}$', fontsize=16)
ax4.set_ylabel(r'$-\sigma_{differential}$ (GPa)', fontsize=16)
ax4.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax4.set_ylim(0.,0.045)

ax5.grid()
#ax3.set_yscale("log")
#ax4.set_ylim(0.,0.1)

ax5.set_xlabel('time', fontsize=16)
ax5.set_ylabel(r'$-\sigma_{reaction and box}$ (GPa)', fontsize=16)
ax5.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
ax5.grid()
#ax3.set_yscale("log")
#ax5.set_ylim(0.,0.03)


#ax6.set_xlabel('ignore plot for now', fontsize=16)
#ax6.set_ylabel(r'ignore plot for now', fontsize=16)
#ax6.legend(bbox_to_anchor=(1.04,1), loc="upper left",fontsize='medium')
#ax6.grid()
##ax3.set_yscale("log")
#ax6.set_ylim(0.,0.1)

fig.tight_layout()
#plt.show()
fig.savefig("geomechanicsTriaxialCompression_TP1_TP2_TP3.png", bbox_inches="tight")

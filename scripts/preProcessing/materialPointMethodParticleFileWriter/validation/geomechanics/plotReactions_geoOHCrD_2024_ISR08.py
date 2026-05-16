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
runLocation='/p/lustre1/malenda1/geosxRuns/'
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

# useDamageAsInitialStrength must be enabled for the painted damage field to be used as strengthScale 
useDamageAsInitialStrength=0               
weibullSeed=1
weibullModulus=6.0
weibullVolume=1.0
flawSize=0.2

sampleHeight = 1. # mm
sampleWidth = 1.  # mm

files=[
'geoOHCrD_2024_ISR08_V5'
# 'stress_control_geoCreepHYD_B-0p08_C-5p7e-10',
]  

shorthand='ISR08'


styles=[
'-',
'--',
'-.',
'-.',
'-.'
]

stressTable=np.array([
[0.,0.,0.,0.],
[0.000000000000000000e+00,-0.000000000000000000e+00,-0.000000000000000000e+00, -0.000000000000000000e+00],
[1.351351351351351426e-02,-8.004315999999999168e-03,-8.004315999999999168e-03, -8.004315999999999168e-03],
[9.864864864864865135e-01,-7.993767000000000580e-03,-7.993767000000000580e-03, -7.993767000000000580e-03],
[1.000000000000000000e+00,1.86157999999999997e-06,1.86157999999999997e-06, 1.86157999999999997e-06]
])

temperatureTable=np.array([
[(0.00e+00), 2.03e+01+273.15],
[(1.51e-01), 2.03e+01+273.15],
[(1.64e-01), 4.08e+01+273.15],
[(3.79e-01), 4.10e+01+273.15],
[(3.92e-01), 6.37e+01+273.15],
[(5.44e-01), 6.40e+01+273.15],
[(5.56e-01), 8.07e+01+273.15],
[(7.72e-01), 8.05e+01+273.15],
[(7.84e-01), 1.03e+02+273.15],
[(9.36e-01), 1.00e+02+273.15],
[(9.49e-01), 4.65e+01+273.15],
[(9.62e-01), 2.63e+01+273.15],
[(9.87e-01), 2.10e+01+273.15]
])


driverTime=stressTable[:,0]
driverPressure=-1.*stressTable[:,1]

driverTemperature=temperatureTable[:,1]

#initialize plots
##plot strain versus time data from the experiment
timeData=[]
strainData=[]

with open ('strain_vs_time_ISR08.csv') as csvfile:
	for line in csvfile:
		line = line.strip();
		values = line.split(',');
		
		#for some reason, there are other characters in the first column. Strip them.
		values[0]=re.sub(r'[^0-9.]','',values[0]);

		values[0]=float(values[0]);
		values[1]=float(values[1]);
		timeData.append(values[0]);
		strainData.append(values[1]);

timeDataNorm=[i/timeData[-1] for i in timeData];

# #reduce the length of time and strain
# Dataindex=[i for i in range(25844) if i%100==0]

# timeDataNorm=[timeDataNorm[i] for i in Dataindex]
# strainData=[strainData[i] for i in Dataindex]
##########################################################################


#initialize plots
evenly_spaced_interval = np.linspace(0, 1, len(files))
colors = [cm.rainbow(x) for x in evenly_spaced_interval]

#fig2, (ax4,ax5,ax6) = plt.subplots(3, 1,figsize=(16,12))

for i,file in enumerate(files):
	fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1,figsize=(10, 40))

	print(file)
	reactionFile=runLocation+file+'/boxAverageHistory.csv'
	#Time,Sxx,Syy,Szz,Sxy,Syz,Sxz,Density,Internal Energy,Damage,Porosity,F00,F11,F22
   #  Time, Sxx, Syy, Szz, Syz, Sxz, Sxy, Density, Damage, Internal Energy, Kinetic Energy, epxx, epyy, epzz, epyz, epxz, epxy, volume, temperature, F00, F11, F22
	data = np.genfromtxt(reactionFile, delimiter=',')
	time = data[:,0]
	boxtime=time
	sxx = data[:,1]
	syy = data[:,2]
	szz = data[:,3]
 
	syz = data[:,4]
	sxz = data[:,5]
	sxy = data[:,6]

	rho = data[:,7]
	damage = data[:,8]
	e = data[:,9]
	temperature = data [:,18]
	evp = (data[:,11]+data[:,12]+data[:,13])


	reactionFile=runLocation+file+'/reactionHistory.csv'
	print('this is the reaction file' + reactionFile)
	data = np.genfromtxt(reactionFile, delimiter=',')
	#poro = data[1:,10]
	F00 = data[:,1]
	F11 = data[:,2]
	F22 = data[:,3]
	stopTime = data[-1,0]
	print('this is the file of interest: '+ str(reactionFile))
	print('this is the stopTime: '+ str(stopTime))
	#scaling the time associated with the temperature driver
	TempdriverTime=(temperatureTable[:,0]/temperatureTable[-1,0])*stopTime
	#print('tempdrivertime' + str(TempdriverTime))




	# this hids all non-monotic time entries, so restart data files are cleaned:
	maxt = 0.0
	mask = np.ones(len(time), dtype=bool)
	for ii,t in enumerate(time):
		if (t<=maxt):
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
	damage = damage[mask,...]
	#poro = poro[mask,...]
	F00 = F00[mask,...]
	F11 = F11[mask,...]
	F22 = F22[mask,...]

	print("F00=",F00[0:10])
	print("F11=",F11[0:10])
	print("F22=",F22[0:10])

	p=(-1.0/3.0)*(sxx+syy+szz)
	print ('this is the pressure data: '+str(p))
	ev = -np.log(F00*F11*F22)

	sys.path.append(runLocation+file+"/")
	job = importlib.import_module('pfw_input_'+file)
	tstop = stopTime


	##########################################################################
	#Modeled Strain vs Time Data Interpolation#
	print ('this is the  modeled time class' + str(type(time)))

	time[0]=0.000
	ev[0]=0.000
	print('time=',time[0:10],'...',time[-1],'], len = ',len(time))
	print('ev=',ev[0:10],'...',ev[-1],'], len = ',len(ev))

	f = interpolate.interp1d(list(time),list(ev))
	print('this is time'+str(time))

	time_Minterpolate=np.arange(time[1],int(tstop),10)
	#Minterpolate means data interplated from the modeled data
	#here, set up array of modeled time data in steps of 10
	ev_Minterpolate=f(time_Minterpolate)

	strainData[0]=0.000
	timeData = tstop*(np.array(timeDataNorm)-timeDataNorm[0])
	print('timeData=',timeData[0:10],'...',timeData[-1],'], len = ',len(timeData))
	print('strainData=',strainData[0:10],'...',strainData[-1],'], len = ',len(strainData))


	g = interpolate.interp1d(list(timeData),list(strainData))
	ev_Dinterpolate=g(time_Minterpolate)

	error=np.sqrt(np.sum(np.array(ev_Minterpolate-ev_Dinterpolate)**2))

	# ##########################################################################
 
 
 
	
	# #Experimental Strain vs Time Data Interpolation#
	# timeDataNorm=np.array(timeDataNorm)
	# #timeDataNorm is the normalized time data from the experiment
	# #here, we are just formatting this data into an array.
	# print('this is time'+str(timeDataNorm))

	# f = interpolate.interp1d(list(tstop*timeDataNorm),list(strainData))
	# #transform the normalized time to regulare time and interpolate experimental strain versus time data

	# strainData_Einterpolate=f(time_Minterpolate)
	# #Einterpolate mean data interplated from the experimental data
	# #interpolate the experimental strain data across the array of modeled time data



	# ##########################################################################
	# #Error Calculation#
	# error =[]
	# for i in range (len(ev_Minterpolate)):
	# 	#for each of the data points in the interpolated time...
		

	# 	errorpt=(abs(ev_Minterpolate[i]-strainData_Einterpolate[i]))**2
	# 	#find the difference between the interpolated modeled and experimental strain data
	# 	#raise to the second power
	# 	error.append(errorpt)
	# error='%.3f' % (sum(error))**0.5
	# #square the sum of all errors at each point.\

	print('this is the error ' + str(error))
	####################################################


		# Compute the equilibrium volumetric plastic strain.
	A = job.creepA
	B = job.creepB
	C = job.creepC / job.timeScale
	D = job.creepD
	E = job.creepE
	F = job.creepF
	p3 = job.p3
	phii = 1. - np.exp(-p3)
	phie = A*np.exp(-driverPressure/B)    # equilibrium unloaded porosity
	print ("this is the equilibrium Phi" + str(phie))
	evpe = np.log( (phii-1.) / (phie-1.) )# equilibirum vol plastic strain
	# estimate the elastic vol strain.
	b0 = job.b0
	b1 = job.b1
	b2 = job.b2
	bulk = b0 + 0.5*b1*np.exp(b2/(-3.*(driverPressure+0.00001))) 
	eve = driverPressure/bulk

	Kd = job.stressControlKd
	Ki = job.stressControlKi
	Kp = job.stressControlKp
 
	Q = job.Q

	######################################################### plotting #########################################################


	ax1.plot(time,p,linestyle='-',color='r',linewidth=2,label='modeled pressure')
	ax1.plot(tstop*driverTime,driverPressure,'--k', label = 'stress table data')
	ax1.legend()


	ax2.plot(time,ev,linestyle='-',color='b',linewidth=2,label='modeled comp. vol. strain')
	ax2.plot(tstop*driverTime,-evpe+eve,linestyle='--',color='c',linewidth=2,label='modeled equilibrium. vol. strain')
	ax2.plot(time_Minterpolate,ev_Minterpolate,'--y',linewidth=2,label='interpolated modeled data')
	ax2.plot(time_Minterpolate,ev_Dinterpolate,'--k',linewidth=2,label='interpolated experimental data')
	ax2.plot(time,ev,linestyle='-',color='k',linewidth=4,alpha=0.5 ,label= 'modeled evp actual versus strain')
	ax2.legend(loc = "lower right")
 

	ax3.plot(ev,p,linestyle='-',color='k',linewidth=2, label= 'modeled pressure versus strain')
	ax3.legend()
 
 
	ax4.plot(ev,temperature,linestyle='-',color='k',linewidth=2, label= 'modeled temperature versus strain')
	ax4.legend()

 
	ax5.plot(TempdriverTime,driverTemperature,'k', label = 'temp table data')
	ax5.plot(boxtime,temperature,'--k', label = 'box temp data')
	ax5.legend()



	######################################################### set labels #########################################################

	ax1.set_xlabel('Stop Time')
	ax1.set_ylabel('Pressure (GPa)')
	ax1.grid()

	ax2.set_xlabel('Stop Time')
	ax2.set_xlim(0.0,tstop)
	ax2.set_ylabel('Compressive Vol. Strain (-)')
	#ax2.set_ylim(0.00,0.06)
	ax2.grid()

	ax3.set_xlabel('Compressive Vol. Strain (-)')
	ax3.set_ylabel('Pressure (GPa)')
	ax3.grid()
 
	ax4.set_xlabel('Compressive Vol Strain (-)')
	ax4.set_ylabel('Temperature (Kelvin)')
	#ax4.set_ylim(-0.2,0.2)
	ax4.grid()
 
 	# von Mises vs density
	ax5.set_xlabel('Stop Time')
	ax5.set_ylabel('Temperature (Kelvin)')
	ax5.grid()


	phi_slope = '-0.0003'
	Q_val = Q
	phi_fit = 'linear'

	fig.suptitle('sample= '+shorthand+'_OHCrD; tstop=' +str(tstop) + '; kp='+str(Kp)+ '; ki='+str(Ki)+'; kd='+str(Kd)+'; A='+str(A)+'; B='+str(B)+ '; C='+str(C)+'D='+str(D)+'E='+str(E)+'F='+str(F)+'; p3=' +str(p3)+'; error=' +str(round(error,3))+'phi_slope'+str(phi_slope)+'; Q:'+str(Q_val)+'phi_fit'+str(phi_fit))
	#fig.suptitle('tstop = 1000; kp = 1.00, ki = 1.00, kd =0.00, A= ' +str(A))
	fig.tight_layout()
	fig.savefig('OHCrD_'+shorthand+'_2024_fit_A='+str(A)+'_B='+str(B)+'_C='+str(C)+'_D='+str(D)+'_E='+str(E)+'_F='+str(F)+'_p3=' +str(p3)+'phi_slope'+str(phi_slope)+'; Q:'+str(Q_val)+'_phi_fit'+str(phi_fit)+'_polynomial.png', bbox_inches="tight")

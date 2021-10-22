import sys
import matplotlib
import numpy as np

import math
import matplotlib.pyplot as plt
import HydrofactureSolutions




if len(sys.argv) < 9:
    sys.exit('Usage: %s prefix mu E\' q KI r_source aper_cutoff' % sys.argv[0])

mu = float(sys.argv[1])
E  = float(sys.argv[2])
q  = float(sys.argv[3])
KI = float(sys.argv[4])
x_source = float(sys.argv[5])
aper_cutoff = float(sys.argv[6])
P0 = float(sys.argv[7])

numFiles = len(sys.argv) - 8
#print numFiles
prefix = []
t_sim = []
radP = []
radA = []
aper0 = []
press0 = []
Area = []
Radius=[]
Pressure = []
Aperture = []
labels = []
symbols = ['ko','rs','k^','gs','rD','bx']
for i in range(0,numFiles):
    prefix.append(sys.argv[8+i])

    t_sim.append(np.empty([0]))
    aper0.append(np.empty([0]))
    press0.append(np.empty([0]))
    Area.append(np.empty([0]))
    Radius.append(np.empty([0]))
    radP.append(np.empty([0]))
    radA.append(np.empty([0]))
    Pressure.append(np.empty([0]))
    Aperture.append(np.empty([0]))

    t_sim[i], press0[i], aper0[i], Area[i] = np.loadtxt(prefix[i]+'_timehist.txt',skiprows=1, unpack=True)
    
    radP[i], Pressure[i] = np.loadtxt(prefix[i]+'_pressure.ply',skiprows=11, usecols = (0,1), unpack=True)
    radA[i], Aperture[i] = np.loadtxt(prefix[i]+'_aperture.ply',skiprows=11, usecols = (0,1), unpack=True)

    #print Aperture[i]
    cutoffList = []
    for j in range(0,len(radP[i])):
        if Aperture[i][j] < aper_cutoff:
            cutoffList.append(j)

    #print cutoffList
    radP[i] = np.delete( radP[i], cutoffList )
    radA[i] = np.delete( radA[i], cutoffList )
    Aperture[i] = np.delete( Aperture[i], cutoffList )
    Pressure[i] = np.delete( Pressure[i], cutoffList )

    labels.append(prefix[i])

    Radius[i] = (Area[i]/math.pi)**0.5
    
#print Aperture[0]
    
#labels[0] = "GEOS Results $K_{IC}=2.0e6$, $\mu=0.001$"
#labels[1] = "GEOS Results $K_{IC}/10$"
#labels[2] = "GEOS Results $\mu/100$"

endTime = t_sim[0][-1]
time = np.arange(0,endTime+1,0.2)
radTimes = np.array([t_sim[0][-1]])
hfsolns = HydrofactureSolutions.PennySolutions()
solns = hfsolns.Solutions( mu, E, q, KI, time, radTimes, x_source )

rho       = solns[0]
w_tough   = solns[1]
p_tough   = solns[2]
w_viscous = solns[3]
p_viscous = solns[4]
P_Tdom    = solns[5]
R_Tdom    = solns[6]
w0_Tdom   = solns[7]
P_Vdom    = solns[8]
R_Vdom    = solns[9]
w0_Vdom   = solns[10]


N1 = 1
font = {'size'   : 10}

matplotlib.rc('font', **font)

fig1 = plt.figure(figsize=(8, 5))  # a new figure window
fig1.subplots_adjust(wspace=0.4,bottom=0.1,left=0.1,right=0.97,top=0.95)




ax2 = fig1.add_subplot(3, 2, 1)  # specify (nrows, ncols, axnum)
plt.plot(time,R_Tdom, 'k-', label='asymptotic $\mu \Rightarrow 0$  soln' )
plt.plot(time,R_Vdom, 'r-', label='asymptotic $K_{IC} \Rightarrow 0$ soln' )
for i in range(0,numFiles):
    plt.plot(t_sim[i][1::N1],Radius[i][1::N1],symbols[i],fillstyle='none', markersize=4, label=labels[i] )
#plt.xlabel('time (s)')
plt.ylabel('Fracture\n Length (m)', multialignment='center')
#plt.xticks(np.arange(min(t_sim[0]), max(t_sim[0]), 20.0))
plt.xlim([0, max(t_sim[0]) ])
plt.legend( bbox_to_anchor=(2.4, 1.1),prop={'size':10})


ax1 = fig1.add_subplot(3, 2, 3)  # specify (nrows, ncols, axnum)
plt.plot(time,w0_Tdom*1000, 'k-', label='$\mu$ => 0' )
plt.plot(time,w0_Vdom*1000, 'r-', label='$K_{IC}$ => 0' )
for i in range(0,numFiles):
    plt.plot(t_sim[i][1::N1],aper0[i][1::N1]*1000,symbols[i],fillstyle='none', markersize=4, label=labels[i] )
#plt.xlabel('time (s)')
plt.ylabel('Aperture(@L=0.5m) \n (mm)', multialignment='center')
#plt.legend(loc='lower right')
#plt.xticks(np.arange(min(t_sim[0]), max(t_sim[0]), 20.0))
plt.xlim([0, max(t_sim[0]) ])
#plt.ylim([ 0.0  , 2.0 ])

ax1 = fig1.add_subplot(3, 2, 5)  # specify (nrows, ncols, axnum)
plt.plot(time,P_Tdom/1.0e6, 'k-', label='$\mu$ => 0' )
plt.plot(time,P_Vdom/1.0e6, 'r-', label='$K_{Ic}$ => 0' )
for i in range(0,numFiles):
    plt.plot(t_sim[i][1::N1],(press0[i][1::N1]-P0)/1.0e6,symbols[i],fillstyle='none', markersize=4, label=labels[i] )
plt.xlabel('time (s)')
plt.ylabel('Pressure(@L=0.5m) \n (MPa)', multialignment='center')
#plt.xticks(np.arange(min(t_sim[0]), max(t_sim[0]), 20.0))
plt.xlim([0, max(t_sim[0]) ])
#plt.yticks(np.arange(0.3, 1.01, 0.1))
plt.ylim([ 0.0  , 1 ])



#plt.savefig(prefix + '_vsTime.eps')

N1 = 1
fig2 = fig1 #plt.figure()  # a new figure window
ax2 = fig2.add_subplot(3, 2, 4)  # specify (nrows, ncols, axnum)
for i in range(0,numFiles):
    plt.plot(radA[i]/np.amax(radA[i]),Aperture[i]*1000,symbols[i],fillstyle='none', markersize=4,label=labels[i] )
#plt.plot(solns[0],w_tough[0]*1000,  'k-', label='$\mu$ => 0' )
plt.plot(solns[0],w_viscous[0]*1000,'r-', label='$K_{Ic}$ => 0' )
plt.plot(solns[0],w_tough[0]*1000,'k-', label='$\mu$ => 0' )
#plt.xlabel('Radius (m)')
plt.ylabel('Aperture(t=200s) \n (mm)', multialignment='center')
#plt.legend(loc='upper right')
plt.xlim([0, 1 ])
#plt.ylim([ 0.0  , 2.0 ])
  
#print Pressure[i][1::N1]
#print np.append(Pressure[i][1::N1],Pressure[i][-1])
ax2 = fig2.add_subplot(3, 2, 6)  # specify (nrows, ncols, axnum)
for i in range(0,numFiles):
    plt.plot(radP[i]/np.amax(radP[i]),(Pressure[i]-P0)/1.0e6,symbols[i],fillstyle='none', markersize=4, label=labels[i] )
#plt.plot(solns[0],p_tough[0]/1.0e6,  'k-', label='$\mu$ => 0' )
plt.plot(solns[0],p_viscous[0]/1.0e6,'r-', label='$K_{Ic}$ => 0' )
plt.plot(solns[0],p_tough[0]/1.0e6,'k-', label='$\mu$ => 0' )
plt.xlabel('Length Coordinate') 
plt.ylabel('Pressure(t=200s) \n (MPa)', multialignment='center')
plt.xlim([0, 1 ])
#plt.yticks(np.arange(-0.2, 0.501, 0.1))
#plt.ylim([ -0.3 , 0.31 ])
#  
# plt.savefig('kgd2.eps',dpi=1200)
plt.show()

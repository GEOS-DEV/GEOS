import sys
import matplotlib
import numpy as np

import math
import matplotlib.pyplot as plt


PI = 3.1415

def ReadPlotData( filename, headerlines, x , y ):
    f = open(filename, 'r')

    lines = f.readlines()

    for i in range(0,headerlines):
        lines.pop(0)

#    print lines[0]
    for line in lines:
        vals = line.split(' ')
        x.append(float(vals[0]))
        y.append(float(vals[1]))



def CalculateTheoreticalSolns( mu, E, q, KI,time ):
    # rho is the dimensionless radius with a range (0,1)
    rho = np.arange(999)
    rho = ( rho + 1 ) / 1000.0

    # Omega_k0 is the dimensionless crack opening for the zero-toughness solution k->0
    # Omega_m0 is the dimensionless crack opening for the zero-viscoisty solution m->0 (m is short for mu)
    # PI_m0 is the dimensinoless net pressure  for 

    Omega_k0 = np.empty_like(rho)
    Omega_m0 = np.empty_like(rho)
    PI_m0 = np.empty_like(rho)
    w_tough = np.empty_like(rho)
    p_tough = np.empty_like(rho)

    B1 = 9.269e-2
    A1 = 3.581e-1
    C1 = 6.846e-1
    C2 = 7.098e-2


    tf = time[-1]
    em = pow( 12*mu / E / tf, 1.0/3.0 )
    Lm = pow( E * pow(q,3) * pow(tf,4) / ( 12*mu ), 1.0/9.0 )

    KIp = 4.0 * math.sqrt( 2.0/PI ) * KI
    ek = pow( pow(KIp/E,6) / q / tf, 0.2)
    Lk = pow( q * E * tf / KIp, 0.4)




    for i in range(0,len(rho)):
        ot1 = ( math.sqrt(70.0)/3.0*C1 + 4.0*math.sqrt(5.0)/9.0 * C2 * ( 13*rho[i]-6 ) ) * pow( 1.0 - rho[i], 2.0/3.0 )
        ot2 = B1 * 8.0 / PI *( math.sqrt( 1.0 - rho[i] ) - rho[i] * math.acos( rho[i] ) )
        Omega_m0[i] = ot1 + ot2
        PI_m0[i] = A1 * ( 2.479 - 2.0 / 3.0 / pow( 1.0 - rho[i], 1.0/3.0 ) ) - B1 * ( math.log( rho[i] / 2.0 ) + 1 )
        Omega_k0[i] = pow( 3.0 / 8.0 / PI, 0.2 ) * math.sqrt( 1 - pow( rho[i],2 ) )

    if 0:
        fig = plt.figure()
        plt.plot(rho,Omega_m0 )
        plt.plot(rho,Omega_k0 )
        fig = plt.figure()
        plt.plot(rho,PI_m0 )
    

    Omega_m0 = Omega_m0 * 0.6955

    w_viscous = em * Lm * Omega_m0
    p_viscous = em * E * PI_m0
    w_tough = ek * Lk * Omega_k0
    p_tough = ek * E * 0.3004 * np.ones_like(rho)






    P_Tdom  = np.empty_like( time )
    R_Tdom  = np.empty_like( time )
    w0_Tdom = np.empty_like( time )
    P_Vdom  = np.empty_like( time )
    R_Vdom  = np.empty_like( time )
    w0_Vdom = np.empty_like( time )

    for i in range(0,len(time)):
        t = time[i]
        if t==0.0:
            t=1.0e-3


        em = pow( 12*mu / E / t, 1.0/3.0 )
        Lm = pow( E * pow(q,3) * pow(t,4) / ( 12*mu ), 1.0/9.0 )

        KIp = 4.0 * math.sqrt( 2.0/PI ) * KI
        ek = pow( pow(KIp/E,6) / q / t, 0.2)
        Lk = pow( q * E * t / KIp, 0.4)


        P_Tdom[i] = ek * E * 0.3004
        R_Tdom[i] = Lk * 0.8546
        w0_Tdom[i] = ek * Lk * Omega_k0[0]

        P_Vdom[i]     = em * E * PI_m0[0]
        R_Vdom[i]     = Lm * 0.6955
        w0_Vdom[i] = em * Lm * Omega_m0[0]

    return [ rho, w_tough, p_tough, w_viscous, p_viscous, P_Tdom, R_Tdom, w0_Tdom, P_Vdom, R_Vdom, w0_Vdom ]















if len(sys.argv) < 6:
    sys.exit('Usage: %s prefix mu E\' q KI' % sys.argv[0])

prefix = sys.argv[1]
mu = float(sys.argv[2])
E  = float(sys.argv[3])
q  = float(sys.argv[4])
KI = float(sys.argv[5])

#mu = 0.0005
#E = 3.2e10
#q = 0.02
#KI = 5.47e5
#KI = 4.236e6

queriesDone = 0

t_sim, aper0 = np.loadtxt('CartMesh/'+prefix+'_aper0.curve',skiprows=1, unpack=True)
t_sim, Area = np.loadtxt('CartMesh/'+prefix+'_area.curve',skiprows=1, unpack=True)
radP, Pressure = np.loadtxt('CartMesh/'+prefix+'_pressure.ply',skiprows=11, usecols = (0,1), unpack=True)
radA, Aperture = np.loadtxt('CartMesh/'+prefix+'_aperture.ply',skiprows=11, usecols = (0,1), unpack=True)

t_sim2, aper02 = np.loadtxt('RadialMesh/'+prefix+'_aper0.curve',skiprows=1, unpack=True)
t_sim2, Area2 = np.loadtxt('RadialMesh/'+prefix+'_area.curve',skiprows=1, unpack=True)
radP2, Pressure2 = np.loadtxt('RadialMesh/'+prefix+'_pressure.ply',skiprows=11, usecols = (0,1), unpack=True)
radA2, Aperture2 = np.loadtxt('RadialMesh/'+prefix+'_aperture.ply',skiprows=11, usecols = (0,1), unpack=True)


Radius = np.sqrt( 4 * Area / PI )
Radius2 = np.sqrt( 4 * Area2 / PI )



solns = CalculateTheoreticalSolns( mu, E, q, KI, t_sim )

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

print rho[0], P_Vdom[0]

N1 = 4
font = {'size'   : 10}

matplotlib.rc('font', **font)

fig1 = plt.figure(figsize=(6, 4))  # a new figure window
fig1.subplots_adjust(wspace=0.4,bottom=0.1,left=0.11,right=0.6,top=0.95)


ax1 = fig1.add_subplot(2, 1, 1)  # specify (nrows, ncols, axnum)
plt.plot(t_sim,w0_Tdom*1000, 'b-', label='$\mu$ = 0' )
plt.plot(t_sim,w0_Vdom*1000, 'r-', label='$K_{Ic}$ = 0' )
plt.plot(t_sim[2::N1],aper0[2::N1]*1000,'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
plt.plot(t_sim2[::N1],aper02[::N1]*1000,'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
#plt.xlabel('time (s)')
plt.ylabel('Aperture0 (mm)')
#plt.legend(loc='lower right')
plt.xticks(np.arange(min(t_sim), max(t_sim), 100.0))
plt.xlim([min(t_sim), max(t_sim) ])


ax2 = fig1.add_subplot(2, 1, 2)  # specify (nrows, ncols, axnum)
plt.plot(t_sim,R_Tdom, 'b-', label='$\mu$ => 0' )
plt.plot(t_sim,R_Vdom, 'r-', label='$K_{Ic}$ => 0' )
plt.plot(t_sim[2::N1],Radius[2::N1],'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
plt.plot(t_sim2[::N1],Radius2[::N1],'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
plt.xlabel('time (s)')
plt.ylabel('Radius (m)')
plt.xticks(np.arange(min(t_sim), max(t_sim), 100.0))
plt.xlim([min(t_sim), max(t_sim) ])

plt.legend(loc='center right', bbox_to_anchor=(1.75, 1.1),prop={'size':10})



plt.savefig(prefix + '_vsTime.eps')
#plt.show()

N1 = 8
fig2 = plt.figure(figsize=(6, 4))  # a new figure window
fig2.subplots_adjust(wspace=0.4,bottom=0.1,left=0.11,right=0.6,top=0.95)

ax2 = fig2.add_subplot(2, 1, 1)  # specify (nrows, ncols, axnum)
plt.plot(solns[0],w_tough*1000,  'b-', label='$\mu$ => 0' )
plt.plot(solns[0],w_viscous*1000,'r-', label='$K_{Ic}$ => 0' )
plt.plot(radA[2::N1]/np.amax(Radius),Aperture[2::N1]*1000,'ks',fillstyle='none', markersize=3,label='GEOS-Cartesian Mesh' )
plt.plot(radA2/np.amax(Radius2),Aperture2*1000,'ko',fillstyle='none', markersize=3,label='GEOS-Radial Mesh' )
#plt.xlabel('Radius (m)')
plt.ylabel('Aperture (mm)')
#plt.legend(loc='upper right')
plt.xlim([0, max(radA/np.amax(Radius)) ])


ax2 = fig2.add_subplot(2, 1, 2)  # specify (nrows, ncols, axnum)
plt.plot(solns[0],p_tough/1.0e6,  'b-', label='$\mu$ = 0' )
plt.plot(solns[0],p_viscous/1.0e6,'r-', label='$K_{Ic}$ = 0' )
plt.plot(radP[2::N1]/np.amax(Radius),Pressure[2::N1]/1.0e6,'ks',fillstyle='none', markersize=3, label='GEOS-Cartesian Mesh' )
plt.plot(radP2[2::N1]/np.amax(Radius2),Pressure2[2::N1]/1.0e6,'ko',fillstyle='none', markersize=3, label='GEOS-Radial Mesh' )
plt.xlabel('Radial Coordinate')
plt.ylabel('Pressure (MPa)')
plt.xlim([0, max(radP/np.amax(Radius)) ])
plt.legend(loc='center right', bbox_to_anchor=(1.75, 1.1),prop={'size':10})



#ylim = plt.ylim()
#pmin = np.amin(Pressure/1.0e6)
#pmax = np.amax(Pressure/1.0e6)
#range = np.ptp(Pressure/1.0e6)
#plt.ylim([ pmin - range/5  , pmax + range/5 ])
#plt.ylim([ -0.1  , 0.4 ])

plt.savefig(prefix + '_vsRadius.eps',dpi=400)


#plt.show()

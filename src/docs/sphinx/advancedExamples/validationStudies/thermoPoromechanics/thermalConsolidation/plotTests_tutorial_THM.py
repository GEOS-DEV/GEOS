import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import h5py
import xml.etree.ElementTree as ElementTree
from mpmath import *
import math
import os


def readGeosTotalStress(depths=(0.0, 4.2, 5.6), biot=1.0):
    """Contrainte TOTALE GEOS aux profondeurs demandees, lue UNIQUEMENT depuis les HDF5
    (sortie TimeHistory de GEOS ; aucune dependance a vtk).
    averageStress = contrainte effective thermo-poro (terme -3aKdT inclus sur cette branche)
    -> contrainte totale = averageStress - biot * pression.
    Necessite stressHistory.hdf5 et pressureHistory.hdf5."""
    if not os.path.exists("stressHistory.hdf5"):
        raise FileNotFoundError(
            "stressHistory.hdf5 introuvable : activez la sortie TimeHistory de averageStress "
            "dans le XML (Task 'stressCollection' + Output 'stressHistoryOutput') puis relancez GEOS.")
    hs = h5py.File("stressHistory.hdf5", 'r')
    ts = np.array(hs.get('averageStress Time')).ravel()
    center = np.array(hs.get('averageStress elementCenter'))   # (nt, ncell, 3)
    S = np.array(hs.get('averageStress'))                      # (nt, ncell, 6)
    hp = h5py.File("pressureHistory.hdf5", 'r')
    P = np.array(hp.get('pressure'))                           # (nt, ncell)
    # coupe au premier reset de temps (redemarrage/restart)
    last = len(ts)
    for j in range(1, len(ts)):
        if ts[j] < 1e-12:
            last = j; break
    n = min(last, len(P))
    ts, S, P = ts[:n], S[:n], P[:n]
    yc = center[0, :, 1]
    cid = {d: int(np.argmin(np.abs(yc - d))) for d in depths}
    out = {}
    for d, i in cid.items():
        out[d] = dict(syy=S[:, i, 1] - biot * P[:, i],
                      sxx=0.5 * (S[:, i, 0] + S[:, i, 2]) - biot * P[:, i])
    return ts, out


class thermalAnalytical:

    def __init__(self, hydromechanicalParameters, xMin, xMax, appliedLoad):
        E = hydromechanicalParameters["YoungModulus"]
        nu = hydromechanicalParameters["PoissonRatio"]
        b = hydromechanicalParameters["biotCoefficient"]
        mu = hydromechanicalParameters["fluidViscosity"]
        cf = hydromechanicalParameters["fluidCompressibility"]
        phi = hydromechanicalParameters["porosity"]
        k = hydromechanicalParameters["permeability"]
        cc = hydromechanicalParameters["consolidationCoefficient"]

        G = E / 2.0/ (1.0 + nu)
        eta = b*(1.-2.*nu)/2./(1.-nu)


        K = E / 3.0 / (1.0 - 2.0 * nu)    # bulk modulus
        Kv = E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))    # uniaxial bulk modulus
        Se = (b - phi) * (1.0 - b) / K + phi * cf    # constrained specific storage

        print( abs((k / mu) * Kv / (Se * Kv + b**2)/cc-1.0) )

        self.characteristicLength = xMax - xMin
        self.appliedLoad = abs(appliedLoad)
        self.loadingEfficiency = b / (Kv * Se + b**2)
        self.consolidationCoefficient = (k / mu) * Kv / (Se * Kv + b**2)
        #self.consolidationCoefficient = cc
        self.initialPressure = 0.0
        self.initialDisplacement = 0.0
        self.term = self.appliedLoad * self.characteristicLength * eta/G

        print( self.term, eta, G )

    def computePressure(self, x, t):
        if t == 0.0:
            return self.initialPressure
        else:
            cc = self.consolidationCoefficient
            L = self.characteristicLength
            load = self.appliedLoad
            p = nsum(
                lambda m: 1 / (2 * m + 1) * exp(-((2 * m + 1)**2) * (math.pi**2) * cc * t / 4 / L / L) * sin(
                    (2 * m + 1) * math.pi * x / 2 / L), [0, inf])
            return load - load* 4 / math.pi * p

    def computeDisplacement(self, x, t):
        if t == 0.0:
            return self.initialDisplacement
        else:
            cc = self.consolidationCoefficient
            L = self.characteristicLength
            term = self.term
            p = nsum(
                lambda m: 8 / ((2 * m + 1)**2) / (math.pi**2) * (exp(-((2 * m + 1)**2) * (math.pi**2) * cc * t / 4 / L / L) -1.0) * cos((2 * m + 1) * math.pi * x / 2 / L), [0, inf])
            return term * p


def getHydromechanicalParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)

    param1 = tree.find('Constitutive/ElasticIsotropic')
    param2 = tree.find('Constitutive/BiotPorosity')
    param3 = tree.find('Constitutive/CompressibleSinglePhaseFluid')
    param4 = tree.find('Constitutive/ConstantPermeability')

    hydromechanicalParameters = dict.fromkeys(["YoungModulus",
                                               "PoissonRatio",
                                               "biotCoefficient",
                                               "fluidViscosity",
                                               "fluidCompressibility",
                                               "porosity",
                                               "permeability",
                                               "skemptonCoefficient",
                                               "poissonRatio",
                                               "undrainedPoissonRatio",
                                               "consolidationCoefficient"])

    hydromechanicalParameters["YoungModulus"] = float(param1.get("defaultYoungModulus"))
    hydromechanicalParameters["PoissonRatio"] = float(param1.get("defaultPoissonRatio"))

    E = hydromechanicalParameters["YoungModulus"]
    nu = hydromechanicalParameters["PoissonRatio"]
    K = E / 3.0/ (1.0 - 2.0 * nu)
    G = E / 2.0/ (1.0 + nu)
    Ks = float(param2.get("grainBulkModulus"))

    hydromechanicalParameters["biotCoefficient"] = 1.0 - K / Ks
    hydromechanicalParameters["porosity"] = float(param2.get("defaultReferencePorosity"))
    hydromechanicalParameters["fluidViscosity"] = float(param3.get("defaultViscosity"))
    hydromechanicalParameters["fluidCompressibility"] = float(param3.get("compressibility"))

    perm = param4.get("permeabilityComponents")
    perm = np.array(perm[1:-1].split(','),float)
    hydromechanicalParameters["permeability"] = perm[0]

    phi = hydromechanicalParameters["porosity"]
    cf = hydromechanicalParameters["fluidCompressibility"]
    bBiot = hydromechanicalParameters["biotCoefficient"]
    kp = hydromechanicalParameters["permeability"]
    mu = hydromechanicalParameters["fluidViscosity"]
    M = 1./(phi*cf + (bBiot - phi)/Ks)
    Ku = K + bBiot**2*M
    B = bBiot*M/Ku
    nuu = (3.*nu + bBiot* B* (1-2.*nu))/(3.-bBiot*B*(1-2.*nu))
    cc = 2.*kp/mu*B**2*G*(1.-nu)*(1.+nuu)**2/9./(1.-nuu)/(nuu-nu)
    hydromechanicalParameters["skemptonCoefficient"] = B
    hydromechanicalParameters["poissonRatio"] = nu
    hydromechanicalParameters["undrainedPoissonRatio"] = nuu
    hydromechanicalParameters["consolidationCoefficient"] = cc 

    return hydromechanicalParameters


def getAppliedTractionFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    param = tree.findall('FieldSpecifications/FieldSpecification')
    load = np.empty(1)
    for elem in param:
        if elem.get("name") == "boundaryPressure" and elem.get("fieldName") == "pressure":
            load = float(elem.get("scale"))
    print(load)

    return load


def getDomainMaxMinXCoordFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath)
    meshElement = tree.find('Mesh/InternalMesh')
    nodeXCoords = meshElement.get("xCoords")
    nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]
    xMin = nodeXCoords[0]
    xMax = nodeXCoords[-1]
    return xMin, xMax


def computePP(x, t, h, phi12, phi11, r11, r22):
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h) * ( phi12 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) + phi11 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) )* sin(((2 * m + 1) * math.pi / 2 / h) * ( h - x) ), [0, inf])
    return 2/h * p

def computeTemp(x, t, h, phi22, phi21, r11, r22):
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h) * ( phi22 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) + phi21 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) )* sin(((2 * m + 1) * math.pi / 2 / h) * ( h - x) ), [0, inf])
    return 2/h * p

def computeDisp(x, t, h, t1, t2, r11, r22, A, M):  
    p = nsum(lambda m: 1 / ((2 * m + 1) * math.pi / 2 / h)**2 * ( t1 * exp(-r11 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) + t2 * exp(-r22 * ((2 * m + 1) * math.pi / 2 / h)**2 * t) )* cos(((2 * m + 1) * math.pi / 2 / h) * ( h - x) ), [0, inf])
    return -A*x+2/h/M * p

def main():
    E = 6.0e3
    nu = 0.4
    Ks = 1.0e27
    Kw = 1.0e50
    a_s = 3.0e-7
    a_w = 0.0
    phi = 0.2
    cs = 1.672e5
    cw = 1.672e2
    rhow = 1.0e3
    rhos = 2.4e3
    Kt = 836

    lamda = E*nu/(1.0+nu)/(1.0-2.0*nu)
    G = E/2.0/(1.0+nu)
    K = E/3.0/(1.0-2.0*nu)
    alpha = 1.0 - K / Ks    
    M = lamda + 2.0*G
    m = (1.0-phi)*cs+phi*cw
    a_m = (1.0-phi)*a_s+phi*a_w
    
    a_m = 9.0e-7
    m = 1.672e5

    beta = a_m*(3.0*lamda+2.0*G)/3.
    #a_p = phi/Kw + (alpha- phi)/Ks
    a_p = 0
    perm = 4.0e-9
    mu = 1.0e-3
    km = perm/mu
    h = 7.0
    theta_0 =0.
    theta_a = 50.
    dT0 = theta_a - theta_0
    F = 1.0

    nuu = 0.499
    M_biot = 2.0*G*(nuu-nu)/(alpha**2*(1.0-2.0*nuu)*(1.0-2.0*nu))
    a_p = 1.0/M_biot
    p0 = alpha*M_biot/(lamda+2.0*G+alpha**2*M_biot)*F
    A = (F-alpha*0-beta*theta_a)/M
    print(p0, A)

    keppa = Kt/m
    c = km*M
    time = np.logspace(-2, 5, 100, endpoint=True)
    T = keppa*time/h**2
    print(a_m, M)
    print(K, alpha, a_p, keppa, M, c/keppa)
    print(K,G,M,A)
   
    g11=1.0/km*(a_p + alpha**2/M)
    g12=1.0/km*(alpha*beta/M-a_m)
    g21=alpha*beta*(theta_0)/(Kt*M)
    g22=(m+beta**2*(theta_0)/M)/Kt
    gs = g11*g22 - g12*g21
    r11=g11/gs
    r12=g12/gs
    r21=g21/gs
    r22=g22/gs
    r1s = r12*(theta_0 - theta_a) - r22*p0
    r2s = r21*p0 - r11*(theta_0 - theta_a)
    rs = r11*r22-r12*r21

    phi12= -(r12*r2s)/rs
    phi11= -(r11*r1s)/rs
    phi22= -(r22*r2s)/rs
    phi21= -(r21*r1s)/rs
    t1 = alpha*phi12 + beta*phi22
    t2 = alpha*phi11 + beta*phi21  

    
    # File paths
    hdf5FilePathPressure = "pressureHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathPressure, 'r')
    timePressure = hf.get('pressure Time')
    timePressure = np.array(timePressure)
    centerPressure = hf.get('pressure elementCenter')
    centerPressure = np.array(centerPressure)
    pressure = hf.get('pressure')
    pressure = np.array(pressure)

    posElement1 = -1
    posElement2 = -1 
    posElement3 = -1
    last = -1
    for j in range(0, centerPressure.shape[1]):
        if centerPressure[0,j,1] >= 0 and posElement1 == -1:
            posElement1 = j
        if centerPressure[0,j,1] >= 4.2 and posElement2 == -1:
            posElement2 = j
        if centerPressure[0,j,1] >= 5.6 and posElement3 == -1:
            posElement3 = j


    # File paths
    hdf5FilePathTemperature = "temperatureHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathTemperature, 'r')
    timeTemperature = hf.get('temperature Time')
    timeTemperature = np.array(timeTemperature)
    centerTemperature = hf.get('temperature elementCenter')
    centerTemperature = np.array(centerTemperature)
    temperature = hf.get('temperature')
    temperature = np.array(temperature)
    
    posElement1 = -1
    posElement2 = -1 
    posElement3 = -1
    last = -1
    
    for j in range(0, centerTemperature.shape[1]):
        if centerTemperature[0,j,1] >= 0 and posElement1 == -1:
            posElement1 = j
        if centerTemperature[0,j,1] >= 4.2 and posElement2 == -1:
            posElement2 = j
        if centerTemperature[0,j,1] >= 5.6 and posElement3 == -1:
            posElement3 = j
 
    for j in range(0, timeTemperature.shape[0]):            
        if j > 0 and timeTemperature[j] < 1e-12 and last == -1:     
            last = j




    #Visulization
    N1 = 4
    fsize = 32
    msize = 8
    lw = 4
    mew = 2
    malpha = 0.6
    lalpha = 0.8

    fig, ax = plt.subplots(2, 2, figsize=(32, 18))
    cmap = plt.get_cmap("tab10")

    pressure_analytical = np.empty(len(time))
    temperature_analytical = np.empty(len(time))
    displacement_analytical = np.empty(len(time))
    x_analytical = [0.0, 0.6, 0.8]
    x_analytical2 = [0.2, 0.6, 1.0]
    iplt = -1
    for xCell in x_analytical:
        iplt += 1
        i = 0
        for k in range(0, len(time)):
            pressure_analytical[i] = computePP(xCell*h, time[i], h, phi12, phi11, r11, r22) 
            temperature_analytical[i] = theta_a + computeTemp(xCell*h, time[i], h, phi22, phi21, r11, r22) 
            displacement_analytical[i] = computeDisp(x_analytical2[iplt]*h, time[i], h, t1, t2, r11, r22, A, M)          
            i += 1            
        #ax[0,0].plot(pressure[k, :]/1.0e6, x[k, :, 0], 'o', color=cmap(iplt), markersize=msize, alpha=malpha, mec=cmap(iplt), fillstyle='none', mew=mew, label='GEOS: t ='+ str(t) + ' s')
        test= np.linspace(0, h, 10, endpoint=True)
        for tt in range(0, len(test)):
            print(computeDisp(test[tt], 0, h, t1, t2, r11, r22, A, M))
        ax[0,0].semilogx(time, pressure_analytical, color=cmap(iplt), lw=lw, alpha=lalpha, label='Analytical: z/h =' + str(xCell))
        #ax[0,1].plot(displacement[k, :, 0], xl_node[k, :, 0], 'o', color=cmap(iplt), markersize=msize, alpha=malpha, mec=cmap(iplt), fillstyle='none', mew=mew, label='GEOS: t ='+ str(t) + ' s')
        ax[0,1].semilogx(time, temperature_analytical+273.15, color=cmap(iplt), lw=lw, alpha=lalpha, label='Analytical: z/h =' + str(xCell))
        ax[1,0].semilogx(time, -displacement_analytical, color=cmap(iplt), lw=lw, alpha=lalpha, label='Analytical: z =' + str(round(x_analytical2[iplt]*h,2)))


    ax[0,0].plot( timePressure[0:last:N1],
              pressure[0:last:N1,posElement1],        
              'o', color=cmap(0), markersize=msize, alpha=malpha, mec=cmap(0), fillstyle='none', mew=mew,
              label='GEOS: z = 0.0 m')
    ax[0,0].plot( timePressure[0:last:N1],
              pressure[0:last:N1,posElement2],        
              'o', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew,
              label='GEOS: z = 4.2 m')
    ax[0,0].plot( timePressure[0:last:N1],
              pressure[0:last:N1,posElement3],        
              'o', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew,
              label='GEOS: z = 5.6 m')

    ax[0,0].set_xlabel('Time [s]', size=fsize, weight="bold")
    ax[0,0].set_ylabel('Pore Pressure [pa]', size=fsize, weight="bold")
    ax[0,0].legend(loc='lower left', fontsize=fsize * 0.6)
    ax[0,0].grid(True)
    ax[0,0].xaxis.set_tick_params(labelsize=fsize)
    ax[0,0].yaxis.set_tick_params(labelsize=fsize) 
    #ax[0,0].invert_yaxis() 


    ax[0,1].plot( timeTemperature[0:last:N1],
              temperature[0:last:N1,posElement1],        
              'o', color=cmap(0), markersize=msize, alpha=malpha, mec=cmap(0), fillstyle='none', mew=mew,
              label='GEOS: z = 0.0 m')
    ax[0,1].plot( timeTemperature[0:last:N1],
              temperature[0:last:N1,posElement2],        
              'o', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew,
              label='GEOS: z = 4.2 m')
    ax[0,1].plot( timeTemperature[0:last:N1],
              temperature[0:last:N1,posElement3],        
              'o', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew,
              label='GEOS: z = 5.6 m')

    ax[0,1].set_xlabel('Time [s]', size=fsize, weight="bold")
    ax[0,1].set_ylabel('Temperature [K]', size=fsize, weight="bold")
    ax[0,1].legend(loc='upper left', fontsize=fsize * 0.6)
    ax[0,1].grid(True)
    ax[0,1].xaxis.set_tick_params(labelsize=fsize)
    ax[0,1].yaxis.set_tick_params(labelsize=fsize) 


    # File paths
    hdf5FilePathDisplacement = "displacementHistory.hdf5"

    
    # Read simulation output from HDF5 file
    hf = h5py.File(hdf5FilePathDisplacement, 'r')
    timeDisplacement = hf.get('totalDisplacement Time')
    timeDisplacement = np.array(timeDisplacement)
    centerDisplacement = hf.get('totalDisplacement ReferencePosition')
    centerDisplacement = np.array(centerDisplacement)
    displacement = hf.get('totalDisplacement')
    displacement = np.array(displacement)

    posVertex1 = -1
    posVertex2 = -1 
    posVertex3 = -1
    last = -1
    
    for j in range(0, 284):            
        if centerDisplacement[0,j,1] >= 1.4 and posVertex1 == -1:
            posVertex1 = j
        if centerDisplacement[0,j,1] >= 4.2 and posVertex2 == -1:
            posVertex2 = j
        if centerDisplacement[0,j,1] >= 7 and posVertex3 == -1:
            posVertex3 = j

    for j in range(0, timeDisplacement.shape[0]):            
        if j > 0 and timeDisplacement[j] < 1e-12 and last == -1:     
            last = j

    ax[1,0].plot( timeDisplacement[0:last:N1],
              -displacement[0:last:N1,posVertex1,1],        
              'o', color=cmap(0), markersize=msize, alpha=malpha, mec=cmap(0), fillstyle='none', mew=mew,
              label='GEOS: z = 1.4 m')
    ax[1,0].plot( timeDisplacement[0:last:N1],
              -displacement[0:last:N1,posVertex2,1],
              'o', color=cmap(1), markersize=msize, alpha=malpha, mec=cmap(1), fillstyle='none', mew=mew,
              label='GEOS: z = 4.2 m')
    ax[1,0].plot( timeDisplacement[0:last:N1],
              -displacement[0:last:N1,posVertex3,1],        
              'o', color=cmap(2), markersize=msize, alpha=malpha, mec=cmap(2), fillstyle='none', mew=mew,
              label='GEOS: z = 7.0 m')

    ax[1,0].set_xlabel('Time [s]', size=fsize, weight="bold")
    ax[1,0].set_ylabel('Displacement [m]', size=fsize, weight="bold")
    ax[1,0].legend(loc='upper left', fontsize=fsize * 0.6)
    ax[1,0].grid(True)
    ax[1,0].xaxis.set_tick_params(labelsize=fsize)
    ax[1,0].yaxis.set_tick_params(labelsize=fsize)

    
    # ------------------------------------------------------------------
    # Contrainte TOTALE (bas-droite, ax[1,1])
    #   deformation uniaxiale, equilibre : sigma_yy^tot = -F (constant)
    #   sigma_xx^tot = sigma_zz^tot = -(lamda/M) F - (2G/M)(beta*theta + alpha*p)
    #   avec beta = a_m (3 lamda + 2G)/3 = 3 K a_s (contrainte thermique),
    #        alpha = coefficient de Biot, M = lamda + 2G (module oedometrique).
    #   GEOS : sigma^tot = averageStress - biot * pression (lu depuis stressHistory.hdf5).
    # ------------------------------------------------------------------
    zdepths = [0.0, 4.2, 5.6]
    ts_geos, geosStress = readGeosTotalStress(depths=tuple(zdepths))

    sigyy_analytical = np.empty(len(time))
    sigxx_analytical = np.empty(len(time))
    iplt = -1
    for xCell in x_analytical:
        iplt += 1
        z = xCell * h
        for i in range(0, len(time)):
            p_an = float(computePP(z, time[i], h, phi12, phi11, r11, r22))
            th_an = float(theta_a + computeTemp(z, time[i], h, phi22, phi21, r11, r22))
            sigxx_analytical[i] = -(lamda / M) * F - (2.0 * G / M) * (beta * th_an + alpha * p_an)
            sigyy_analytical[i] = -F
        ax[1, 1].semilogx(time, sigxx_analytical, color=cmap(iplt), lw=lw, alpha=lalpha,
                          label='Analytical $\\sigma_{xx}$: z =' + str(round(z, 2)))
        g = geosStress[zdepths[iplt]]
        ax[1, 1].plot(ts_geos, g['sxx'], 'o', color=cmap(iplt), markersize=msize, alpha=malpha,
                      mec=cmap(iplt), fillstyle='none', mew=mew,
                      label='GEOS $\\sigma_{xx}$: z =' + str(round(z, 2)))

    # contrainte verticale totale (= -F pour toutes les profondeurs)
    ax[1, 1].semilogx(time, -F * np.ones(len(time)), 'k--', lw=lw * 0.6, alpha=0.7,
                      label='$\\sigma_{yy}$ total = -F (all z)')
    ax[1, 1].plot(ts_geos, geosStress[0.0]['syy'], 'k+', markersize=msize, mew=mew, alpha=malpha,
                  label='GEOS $\\sigma_{yy}$ total')

    ax[1, 1].set_xlabel('Time [s]', size=fsize, weight="bold")
    ax[1, 1].set_ylabel('Total stress [Pa]', size=fsize, weight="bold")
    ax[1, 1].legend(loc='upper left', fontsize=fsize * 0.5)
    ax[1, 1].grid(True)
    ax[1, 1].xaxis.set_tick_params(labelsize=fsize)
    ax[1, 1].yaxis.set_tick_params(labelsize=fsize)

    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0.4, hspace=0.4)

    fig.savefig('Verification_tutorial.png')
    print('figure -> Verification_tutorial.png')


if __name__ == "__main__":
    main()

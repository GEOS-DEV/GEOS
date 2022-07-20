import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree

class PKN_viscosityStorageDominated:
    def __init__( self, E, nu, KIC, mu, Q0, t, h ):
        Ep = E / ( 1.0 - nu**2.0 )        
        self.t  = t
        self.Q0 = Q0
        self.mu = mu
        self.Ep = Ep        
        self.h = h

    def analyticalSolution( self ):
        t  = self.t
        Q0 = self.Q0
        mu = self.mu        
        Ep = self.Ep       
        h = self.h

        halfLength = 0.3817 * ( ( Ep * Q0**3.0 * t**4.0 ) / ( mu * h**4.0 ))**( 1.0/5.0 )

        inletAperture = 3.0 * ( ( mu * Q0 * halfLength ) / ( Ep ))**( 1.0/4.0 )        
        
        inletPressure = ( ( 16.0 * mu * Q0 * Ep**3.0 * halfLength ) / ( np.pi * h**4.0 ))**( 1.0/4.0 )
        
        return [ halfLength, inletAperture , inletPressure ]


def getParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse(xmlFilePath + "_benchmark.xml")

    maxTime = float(tree.find('Events').get('maxTime'))

    toughness = float( tree.find('Solvers/SurfaceGenerator').get('rockToughness') )

    param = tree.findall('Geometry/Box')
    found_core = False
    for elem in param:
        if elem.get("name") == "core":
            source = elem.get("xMax")
            source = [float(i) for i in source[1:-1].split(",")]
            height = round(source[1])*2.0
            found_core = True
        if found_core: break
  

    tree = ElementTree.parse(xmlFilePath + "_base.xml")

    elasticParam = tree.find('Constitutive/ElasticIsotropic')
    
    K = float(elasticParam.get('defaultBulkModulus'))
    G = float(elasticParam.get('defaultShearModulus'))
    youngModulus = (9.0 * K * G) / (3.0 * K + G)
    poissonRatio = youngModulus / (2.0 * G) - 1.0    

    viscosity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity') ) 

    fluidDensity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity') )

    injectionRate = -4.0 * float( tree.find('FieldSpecifications/SourceFlux').get('scale') ) / fluidDensity

    return [ maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate, height ]


def main():
    xmlFilePathPrefix = "../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated" 

    tMax, E, nu, KIC, mu, Q0, h = getParametersFromXML( xmlFilePathPrefix )

    t = np.arange( 0.01*tMax, tMax, 0.01*tMax )

    pknFrac = PKN_viscosityStorageDominated ( E, nu, KIC, mu, Q0, t, h )
    halfLength, inletAperture , inletPressure = pknFrac.analyticalSolution()

    # Load GEOSX results
    t_sim, p0_sim, w0_sim, fracArea_sim = np.loadtxt("model-results.txt", skiprows=1, unpack=True)
    halfLength_sim = 2.0*fracArea_sim/h

    # Visulization
    N1 = 1
    fsize = 30
    msize = 12
    lw=8
    mew=2
    malpha = 1.0
    lablelist = ['Analytical Solution', 'GEOSX Result']
    
    fig, ax = plt.subplots(2,2,figsize=(24, 18))
    cmap = plt.get_cmap("tab10")

    ax[0,0].plot(t, halfLength, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0,0].plot(t_sim, halfLength_sim, 's', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[0,0].set_xlim([min(t_sim), max(t_sim)])
    ax[0,0].set_ylim(0, 140)
    ax[0,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,0].set_ylabel(r'Fracture Half Length (m)', size=fsize, weight="bold")
    ax[0,0].legend(loc='lower right',fontsize=fsize*0.7)
    ax[0,0].grid(True)
    ax[0,0].xaxis.set_tick_params(labelsize=fsize)
    ax[0,0].yaxis.set_tick_params(labelsize=fsize)


    ax[0,1].plot(t, inletAperture*1000, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0,1].plot(t_sim, w0_sim*1000, 's', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[0,1].set_xlim([min(t_sim), max(t_sim)])
    ax[0,1].set_ylim(0, 2.0)
    ax[0,1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,1].set_ylabel(r'Fracture Mouth Opening (mm)', size=fsize, weight="bold")
    ax[0,1].legend(loc='lower right',fontsize=fsize*0.7)
    ax[0,1].grid(True)
    ax[0,1].xaxis.set_tick_params(labelsize=fsize)
    ax[0,1].yaxis.set_tick_params(labelsize=fsize)


    ax[1,0].plot(t, inletPressure/1.0e6, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[1,0].plot(t_sim, p0_sim/1.0e6, 's', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[1,0].set_xlim([min(t_sim), max(t_sim)])
    ax[1,0].set_ylim(0, 2.5)
    ax[1,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[1,0].set_ylabel(r'Net Pressure at Well (MPa)', size=fsize, weight="bold")
    ax[1,0].legend(loc='lower right',fontsize=fsize*0.7)
    ax[1,0].grid(True)
    ax[1,0].xaxis.set_tick_params(labelsize=fsize)
    ax[1,0].yaxis.set_tick_params(labelsize=fsize)


    ax[1,1].axis('off')

    plt.show()


if __name__ == "__main__":
    main()

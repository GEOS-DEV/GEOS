import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree

class Penny_toughnessStorageDominated:
    def __init__( self, E, nu, KIC, mu, Q0, t, xSource ):
        Ep = E / ( 1.0 - nu**2.0 )
        Mp = 12.0 * mu
        self.t  = t
        self.Q0 = Q0
        self.mu = mu
        self.Ep = Ep
        self.Mp = Mp
        self.xSource = xSource

    def analyticalSolution( self ):
        t  = self.t
        Q0 = self.Q0
        mu = self.mu        
        Ep = self.Ep
        Mp = self.Mp
        xSource = self.xSource

        def PI_m0(rho):
            A1 = 0.3581
            B = 0.09269
            return A1 * ( 2.479 - 2.0 / 3.0 / pow( 1.0 - rho , 1.0/3.0 )) - B * ( np.log(rho/2.0) + 1.0 )

        fracRadius = 0.6955 * ( ( Ep * Q0**3.0 * t**4.0 ) / ( Mp ))**( 1.0/9.0 )

        inletAperture = 1.1977 * ( ( Mp**2.0 * Q0**3.0 * t ) / ( Ep**2.0 ))**( 1.0/9.0 )        
        
        xd0 = xSource/fracRadius
        inletPressure = PI_m0(xd0) * ( ( Ep**2.0 * Mp ) / ( t ))**( 1.0/3.0 )
        
        return [ fracRadius, inletAperture , inletPressure ]


def getParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse(xmlFilePath + "_benchmark.xml")

    maxTime = float(tree.find('Events').get('maxTime'))

    toughness = float( tree.find('Solvers/SurfaceGenerator').get('rockToughness') )

    param = tree.findall('Geometry/Box')
    found_source = False
    for elem in param:
        if elem.get("name") == "source":
            source = elem.get("xMax")
            source = [float(i) for i in source[1:-1].split(",")]
            x_source = round(source[0])
            found_source = True
        if found_source: break
  

    tree = ElementTree.parse(xmlFilePath + "_base.xml")

    elasticParam = tree.find('Constitutive/ElasticIsotropic')
    
    K = float(elasticParam.get('defaultBulkModulus'))
    G = float(elasticParam.get('defaultShearModulus'))
    youngModulus = (9.0 * K * G) / (3.0 * K + G)
    poissonRatio = youngModulus / (2.0 * G) - 1.0    

    viscosity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity') ) 

    fluidDensity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity') )

    injectionRate = -4.0 * float( tree.find('FieldSpecifications/SourceFlux').get('scale') ) / fluidDensity

    return [ maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate, x_source ]


def main():
    xmlFilePathPrefix = "../../../../../../inputFiles/hydraulicFracturing/pennyShapedViscosityDominated"    

    tMax, E, nu, KIC, mu, Q0, xSource = getParametersFromXML( xmlFilePathPrefix )

    t = np.arange( 0.01*tMax, tMax, 0.01*tMax )

    pennyFrac = Penny_toughnessStorageDominated ( E, nu, KIC, mu, Q0, t, xSource )
    fracRadius, inletAperture , inletPressure = pennyFrac.analyticalSolution()

    # Load GEOSX results
    t_sim, p0_sim, w0_sim, fracArea_sim = np.loadtxt("model-results.txt", skiprows=1, unpack=True)
    Radius_sim = (fracArea_sim/np.pi)**0.5

    # Visulization
    N1 = 1
    fsize = 30
    msize = 12
    lw=8
    mew=2
    malpha = 1.0
    lablelist = ['Asymptotic ($K_{IC}$ => 0, $C_{L}$ => 0 )', 'GEOSX ($K_{IC}$ => 0, $C_{L}$ => 0 )']
    
    fig, ax = plt.subplots(2,2,figsize=(24, 18))
    cmap = plt.get_cmap("tab10")

    ax[0,0].plot(t, fracRadius, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0,0].plot(t_sim, Radius_sim, 'D', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[0,0].set_xlim([min(t_sim), max(t_sim)])
    ax[0,0].set_ylim(0, 80)
    ax[0,0].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,0].set_ylabel(r'Fracture Radius (m)', size=fsize, weight="bold")
    ax[0,0].legend(loc='lower right',fontsize=fsize*0.7)
    ax[0,0].grid(True)
    ax[0,0].xaxis.set_tick_params(labelsize=fsize)
    ax[0,0].yaxis.set_tick_params(labelsize=fsize)


    ax[0,1].plot(t, inletAperture*1000, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[0,1].plot(t_sim, w0_sim*1000, 'D', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[0,1].set_xlim([min(t_sim), max(t_sim)])
    ax[0,1].set_ylim(0, 2.0)
    ax[0,1].set_xlabel(r'Time (s)', size=fsize, weight="bold")
    ax[0,1].set_ylabel(r'Fracture Mouth Opening (mm)', size=fsize, weight="bold")
    ax[0,1].legend(loc='lower right',fontsize=fsize*0.7)
    ax[0,1].grid(True)
    ax[0,1].xaxis.set_tick_params(labelsize=fsize)
    ax[0,1].yaxis.set_tick_params(labelsize=fsize)


    ax[1,0].plot(t, inletPressure/1.0e6, lw=lw, alpha=0.8, color=cmap(0), label=lablelist[0])
    ax[1,0].plot(t_sim, p0_sim/1.0e6, 'D', color=cmap(0), mec = cmap(0), fillstyle='none', markersize=msize, mew=mew, label=lablelist[1], alpha=malpha )
    ax[1,0].set_xlim([min(t_sim), max(t_sim)])
    ax[1,0].set_ylim(0, 1.5)
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

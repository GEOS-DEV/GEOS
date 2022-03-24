import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree

class kgd:
    def __init__(self, E, nu, KIc, mu, Q0, t):
        Ep = E / ( 1.0 - nu**2.0 )

        self.t  = t
        self.Q0 = Q0
        self.mu = mu
        self.Ep = Ep
        self.X  = 256.0 / ( 3.0 * np.pi**2.0 ) * ( KIc**4.0 ) / ( mu * Q0 * Ep**3.0 )


    def analyticalSolution(self):
        t  = self.t
        mu = self.mu
        Q0 = self.Q0
        Ep = self.Ep
        X  = self.X

        halfLength = 0.9324 * X**( -1.0/6.0 ) * ( ( Ep * Q0**3.0 ) / ( 12.0 * mu ) )**( 1.0/6.0 ) * t**( 2.0/3.0 )

        inletAperture = np.sqrt( 0.5 * X**( 0.5 ) * ( 12.0 * mu * Q0 / Ep )**( 0.5 ) * halfLength )

        inletPressure = 0.125 * X**( 0.5 ) * (12.0 * mu * Q0 * Ep)**( 0.5 ) / inletAperture

        return [ halfLength, inletAperture , inletPressure ]

def getParametersFromXML( xmlFilePath ):
    tree = ElementTree.parse(xmlFilePath + "_benchmark.xml")

    maxTime = float(tree.find('Events').get('maxTime'))

    tree = ElementTree.parse(xmlFilePath + "_base.xml")

    elasticParam = tree.find('Constitutive/ElasticIsotropic')

    youngModulus = float(elasticParam.get('defaultYoungModulus'))
    poissonRatio = float(elasticParam.get('defaultPoissonRatio'))

    viscosity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity') )

    toughness = float( tree.find('Solvers/SurfaceGenerator').get('rockToughness') )

    fluidDensity = float( tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity') )

    injectionRate = -2.0 * float( tree.find('FieldSpecifications/SourceFlux').get('scale') ) / fluidDensity

    return [ maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate ]


def main():
    xmlFilePathPrefix = "../../../../../../inputFiles/hydraulicFracturing/kgdToughnessDominated"

    tMax, E, nu, KIc, mu, Q0 = getParametersFromXML( xmlFilePathPrefix )

    t    = np.arange(0.01*tMax,tMax,0.01*tMax)

    model = kgd( E, nu, KIc, mu, Q0, t )
    halfLength, inletAperture , inletPressure = model.analyticalSolution()


    # GEOSX results
    t_geosx, halfLength_geosx = [], []
    for line in open('surfaceArea.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip()==''):
            values = [float(s) for s in line.split()]
            t_geosx.append( values[0] )
            halfLength_geosx.append( values[1] / 2.0 )

    inletAperture_geosx = []
    for line in open('inletAperture.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip()==''):
            values = [float(s) for s in line.split()]
            inletAperture_geosx.append( values[1] * 1e3 )

    inletPressure_geosx = []
    for line in open('inletPressure.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip()==''):
            values = [float(s) for s in line.split()]
            inletPressure_geosx.append( values[1] / 1e6 )


    fig = plt.figure(figsize=[15,10])

    plt.subplot(221)
    plt.plot(t_geosx, halfLength_geosx, 'ko', label='GEOSX result')
    plt.plot(t, halfLength,  'k', linewidth=2, label='Analytic')
    plt.ylabel('Fracture half-length (m)')
    plt.xlabel('Injection time (s)')

    plt.subplot(222)
    plt.plot(t_geosx, inletAperture_geosx, 'ko', label='GEOSX result')
    plt.plot(t, inletAperture * 1e3,  'k', linewidth=2, label='Analytic')
    plt.ylabel('Inlet aperture (mm)')
    plt.xlabel('Injection time (s)')

    plt.subplot(223)
    plt.plot(t_geosx, inletPressure_geosx, 'ko', label='GEOSX result')
    plt.plot(t, inletPressure / 1e6,  'k', linewidth=2, label='Analytic')
    plt.ylabel('Inlet fluid pressure (MPa)')
    plt.xlabel('Injection time (s)')

    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
    
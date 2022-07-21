import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ElementTree
import HydrofractureSolutions


def getParametersFromXML(xmlFilePath):
    tree = ElementTree.parse(xmlFilePath + "_benchmark.xml")

    maxTime = float(tree.find('Events').get('maxTime'))

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

    youngModulus = float(elasticParam.get('defaultYoungModulus'))
    poissonRatio = float(elasticParam.get('defaultPoissonRatio'))

    viscosity = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultViscosity'))

    toughness = float(tree.find('Solvers/SurfaceGenerator').get('rockToughness'))

    fluidDensity = float(tree.find('Constitutive/CompressibleSinglePhaseFluid').get('defaultDensity'))

    injectionRate = -2.0 * float(tree.find('FieldSpecifications/SourceFlux').get('scale')) / fluidDensity

    return [maxTime, youngModulus, poissonRatio, toughness, viscosity, injectionRate, x_source]


def main():
    xmlFilePathPrefix = "../../../../../../../inputFiles/hydraulicFracturing/kgdViscosityDominated"

    tMax, E, nu, KIC, mu, Q0, xSource = getParametersFromXML(xmlFilePathPrefix)
    Ep = E / (1.0 - nu**2.0)

    t = np.arange(0.01 * tMax, tMax, 0.01 * tMax)
    radTimes = np.array([tMax])
    hfsolns = HydrofractureSolutions.KGDSolutions()
    kgdFrac = hfsolns.Solutions(mu, Ep, Q0, KIC, t, radTimes, xSource)
    inletPressure = kgdFrac[8]
    halfLength = kgdFrac[9]
    inletAperture = kgdFrac[10]

    # GEOSX results
    t_geosx, halfLength_geosx = [], []
    for line in open('surfaceArea.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            t_geosx.append(values[0])
            halfLength_geosx.append(values[1] / 2.0)

    inletAperture_geosx = []
    for line in open('inletAperture.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            inletAperture_geosx.append(values[1] * 1e3)

    inletPressure_geosx = []
    for line in open('inletPressure.curve', 'r'):
        if not (line.strip().startswith("#") or line.strip() == ''):
            values = [float(s) for s in line.split()]
            inletPressure_geosx.append(values[1] / 1e6)

    fig = plt.figure(figsize=[15, 10])

    plt.subplot(221)
    plt.plot(t_geosx, halfLength_geosx, 'ko', label='GEOSX result')
    plt.plot(t, halfLength, 'k', linewidth=2, label='Analytical solution')
    plt.ylabel('Fracture half-length (m)')
    plt.xlabel('Injection time (s)')

    plt.subplot(222)
    plt.plot(t_geosx, inletAperture_geosx, 'ko', label='GEOSX result')
    plt.plot(t, inletAperture * 1e3, 'k', linewidth=2, label='Analytical solution')
    plt.ylabel('Inlet aperture (mm)')
    plt.xlabel('Injection time (s)')

    plt.subplot(223)
    plt.plot(t_geosx, inletPressure_geosx, 'ko', label='GEOSX result')
    plt.plot(t, inletPressure / 1e6, 'k', linewidth=2, label='Analytical solution')
    plt.ylabel('Inlet fluid pressure (MPa)')
    plt.xlabel('Injection time (s)')

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()

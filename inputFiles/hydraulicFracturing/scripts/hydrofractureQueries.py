import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin,cos,tan,exp
import hdf5_wrapper
import xml.etree.ElementTree as ElementTree


def getFracHeightFromXML(fracType):
    prefix = "../../../../../../../inputFiles/hydraulicFracturing/"

    tree = ElementTree.parse(prefix + fracType + "_benchmark.xml") 

    if fracType == 'kgdToughnessDominated' or fracType == 'kgdViscosityDominated':
        meshElement = tree.find('Mesh/InternalMesh')
        nodeZCoords = meshElement.get("zCoords")
        nodeZCoords = [float(i) for i in nodeZCoords[1:-1].split(",")]    
        zMin = nodeZCoords[0]
        zMax = nodeZCoords[1]
        height = zMax - zMin

    elif fracType == 'pknViscosityDominated':
        param = tree.findall('Geometry/Box')
        found_core = False
        for elem in param:
            if elem.get("name") == "core":
                source = elem.get("xMax")
                source = [float(i) for i in source[1:-1].split(",")]
                height = round(source[1])
                found_core = True
            if found_core: break

    return height


def main():
    # Load and process GEOSX results
    # File path
    fracType = sys.argv[1] 
    prefix = "../../../../../../../inputFiles/hydraulicFracturing/"   
    hdf5File = prefix + fracType + '_output.hdf5' 

    # Get frac height from XML file
    if fracType == 'kgdToughnessDominated' or fracType == 'kgdViscosityDominated' or fracType == 'pknViscosityDominated':
        fracHeight = getFracHeightFromXML(fracType)

    # Read simulation output from HDF5 file
    # Global Coordinate of Element Center
    hf = hdf5_wrapper.hdf5_wrapper(hdf5File)
    xl = hf['pressure elementCenter']
    xl = np.array(xl)
    xcord = xl[-1,:,0]
    ycord = xl[-1,:,1]
    zcord = xl[-1,:,2]
    tl = hf['pressure Time']
    tl = np.array(tl)
    # Load pressure
    fpre = hf['pressure']  
    fpre = np.array(fpre)
    # Load elementAperture
    aper = hf['elementAperture'] 
    aper = np.array(aper)
    # Load elementArea
    area = hf['elementArea'] 
    area = np.array(area)

    # Find injection location
    ind = np.argmin(ycord)
    r2 = np.zeros([ind])
    for i in range(0,ind):
        r2[i] = xcord[i]**2 + ycord[i]**2 
      
    well_ind = np.argmin(r2)    

    # Compute frac length
    length = np.zeros([len(tl)]) 
    for i in range(0,len(tl)):
        fracArea = 0
        for j in range(0,len(aper[0,:])):
            if aper[i,j]>1.0e-5:
                fracArea += area[i,j]        
        
        if fracType == 'kgdToughnessDominated' or fracType == 'kgdViscosityDominated':
            length[i] = fracArea/fracHeight
        elif fracType == 'pennyShapedToughnessDominated' or fracType == 'pennyShapedViscosityDominated':
            length[i] = (fracArea*4.0/math.pi)**0.5 
        elif fracType == 'pknViscosityDominated':
            length[i] = fracArea/fracHeight


    header   = '      time   pressure   aperture     length'
    timehist = []
    for i in range(0,len(tl)):
        time = tl[i]
        injectionPressure = fpre[i,well_ind]
        injectionAperture = aper[i,well_ind]
        fractureLength = length[i]
        timehist.append([float(time), float(injectionPressure), float(injectionAperture), float(fractureLength)])

    np.savetxt('model-results.txt', timehist, header=header, fmt='%10.4g', comments='')

if __name__ == "__main__":
    main()


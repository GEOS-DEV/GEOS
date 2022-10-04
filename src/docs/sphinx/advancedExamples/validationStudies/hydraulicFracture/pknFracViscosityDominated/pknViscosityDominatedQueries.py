import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from math import sin,cos,tan,exp
import h5py
import xml.etree.ElementTree as ElementTree


def getMeshSettings(xmlFilePath):
    
    tree = ElementTree.parse(xmlFilePath)
    meshElement = tree.find('Mesh/InternalMesh')

    nodeXCoords = meshElement.get("xCoords")
    nodeYCoords = meshElement.get("yCoords")
    nodeZCoords = meshElement.get("zCoords")
    
    nodeXCoords = [float(i) for i in nodeXCoords[1:-1].split(",")]
    nodeYCoords = [float(i) for i in nodeYCoords[1:-1].split(",")]
    nodeZCoords = [float(i) for i in nodeZCoords[1:-1].split(",")]    

    xMin = nodeXCoords[0]
    xMax = nodeXCoords[1]
    yMin = nodeYCoords[0]
    yMax = nodeYCoords[1]
    zMin = nodeZCoords[0]
    zMax = nodeZCoords[1]

    nXElem = meshElement.get("nx")
    nYElem = meshElement.get("ny")
    nZElem = meshElement.get("nz")
    nXElem = [float(i) for i in nXElem[1:-1].split(",")]
    nYElem = [float(i) for i in nYElem[1:-1].split(",")]
    nZElem = [float(i) for i in nZElem[1:-1].split(",")]  
    nx = nXElem[0]
    ny = nYElem[0]
    nz = nZElem[0]
    
    return xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz

def SaveTimeHistory( name, array, header ):
    with open(name, 'w') as filehandle:
        filehandle.write("%s\n" % header )
        for row in array:
            for entry in row:
                filehandle.write("%10.4g" % entry)
            filehandle.write("\n")

def main():
    # Load and process GEOSX results
    # File path
    hdf5File1Path = "../../../../../../../inputFiles/hydraulicFracturing/Penny_zeroToughness_output.hdf5"
    xmlFile1Path = "../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_base.xml"
    xmlFile2Path = "../../../../../../../inputFiles/hydraulicFracturing/pknViscosityDominated_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz = getMeshSettings(xmlFile2Path)
    fracHeight = 20
    dx = abs(xMax-xMin)/nx
    dy = abs(yMax-yMin)/ny

    # Read simulation output from HDF5 file
    # Global Coordinate of Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('pressure elementCenter')
    xl = np.array(xl)
    xcord = xl[0,:,0]
    ycord = xl[0,:,1]
    zcord = xl[0,:,2]
    tl = hf.get('pressure Time')
    tl = np.array(tl)
    # Load pressure
    fpre = hf.get('pressure')  
    fpre = np.array(fpre)
    # Load elementAperture
    aper = hf.get('elementAperture')  
    aper = np.array(aper)
    # Load elementArea
    area = hf.get('elementArea')  
    area = np.array(area)

    # Find injection location
    for i in range(0,len(zcord)):
        if abs(xcord[i]/(dx/2.) - 1) < 0.01 and abs(ycord[i]/(dy/2.) - 1) < 0.01 and abs(zcord[i]) < 0.01:
           well_ind = i

    # Compute frac length
    length = np.zeros([len(tl)]) 
    farea = np.zeros([len(tl)])    
    for i in range(0,len(tl)):
        temp = 0
        for j in range(0,len(aper[0,:])):
            if aper[i,j]>1.0e-5:
               temp += area[i,j]
        
        farea[i] = temp
        length[i] = 2*temp/fracHeight


    header   = [['      time', '  pressure', '  aperture', '    length']]
    timehist = []
    for i in range(0,len(tl)):
        time = tl[i]
        injectionPressure = fpre[i,well_ind]
        injectionAperture = aper[i,well_ind]
        fractureLength = length[i]
        timehist.append([float(time), float(injectionPressure), float(injectionAperture), float(fractureLength)])

    SaveTimeHistory('model-results.txt', timehist, header)

if __name__ == "__main__":
    main()


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
    xMax = nodeXCoords[-1]
    yMin = nodeYCoords[0]
    yMax = nodeYCoords[-1]
    zMin = nodeZCoords[0]
    zMax = nodeZCoords[-1]

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
                filehandle.write("%11.4g" % entry)
            filehandle.write("\n")

def main():
    # Load and process GEOSX results
    # File path
    hdf5File1Path = "../../../../../../../inputFiles/hydraulicFracturing/KGD_validation_output.hdf5"
    xmlFile1Path = "../../../../../../../inputFiles/hydraulicFracturing/kgdValidation_base.xml"
    xmlFile2Path = "../../../../../../../inputFiles/hydraulicFracturing/kgdValidation_benchmark.xml"

    # Read simulation parameters from XML file
    xMin, xMax, yMin, yMax, zMin, zMax, nx, ny, nz = getMeshSettings(xmlFile2Path)
    fracHeight = abs(zMax-zMin)  
    dx = abs(xMax-xMin)/nx  
    dy = abs(yMax-yMin)/ny
    dz = abs(zMax-zMin)/nz
    zmean = (zMax+zMin)/2.     

    # Read simulation output from HDF5 file
    # Global Coordinate of Element Center
    hf = h5py.File(hdf5File1Path, 'r')
    xl = hf.get('pressure elementCenter')
    xl = np.array(xl)
    xcord = xl[-1,:,0]
    ycord = xl[-1,:,1]
    zcord = xl[-1,:,2]
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
    
    # Query simulation results
    xloc_58 = 0.015
    xloc_57 = 0.041
    xloc_lvdt = 0.0285
    wellPressure = np.zeros([len(tl)])
    G58Pressure = np.zeros([len(tl)])
    G57Pressure = np.zeros([len(tl)])
    LVDTAperture = np.zeros([len(tl)])
    length = np.zeros([len(tl)]) 
    farea = np.zeros([len(tl)]) 
    for j in range(0,len(tl)): 
        xcord = xl[j,:,0]
        ycord = xl[j,:,1]
        zcord = xl[j,:,2]
        temp = 0
        for i in range(0,len(aper[0,:])):
            if abs(xcord[i]/(dx/2.) - 1) < 0.01 and abs(ycord[i]) < 0.01 and abs(zcord[i]/(zmean-dz/2.) - 1) < 0.01:
               wellPressure[j] = fpre[j,i]

            if abs(xcord[i]-xloc_58) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i]/(zmean-dz/2.) - 1) < 0.01:
               G58Pressure[j] = fpre[j,i]

            if abs(xcord[i]-xloc_57) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i]/(zmean-dz/2.) - 1) < 0.01:
               G57Pressure[j] = fpre[j,i]

            if abs(xcord[i]-xloc_lvdt) < 0.001 and abs(ycord[i]) < 0.01 and abs(zcord[i]/(zmean-dz/2.) - 1) < 0.01:
               LVDTAperture[j] = aper[j,i]

            if aper[j,i]>1.0e-5:
               temp += area[j,i]
        
        farea[j] = temp        


    # Ouput fracture characteristics
    header   = [['      time', ' wpressure', '58pressure', '57pressure', ' Laperture', '      area']]    
    timehist = []
    for i in range(0,len(tl)):
        timehist.append([float(tl[i]), float(wellPressure[i]), float(G58Pressure[i]), float(G57Pressure[i]), float(LVDTAperture[i]), float(farea[i])])

    SaveTimeHistory('model_results.txt', timehist, header)


if __name__ == "__main__":
    main()


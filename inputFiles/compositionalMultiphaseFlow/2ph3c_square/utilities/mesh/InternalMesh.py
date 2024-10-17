# ------------------------------------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# Copyright (c) 2019-     INRIA project-team Makutu 
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
# ------------------------------------------------------------------------------------------------------------

import numpy as np


class InternalMesh:
    """
    GEOSX Internal Mesh

    Attributes
    ----------
        xml : XML
            XML object containing the information on the mesh
        bounds : list of list
            Real bounds of the mesh [[xmin, xmax],[ymin,ymax],[zmin,zmax]]
        nx : int
            Number of elements in the x direction
        ny : int
            Number of elements in the y direction
        nz : int
            Number of elements in the z direction
        order : int
            Mesh order
        cellBlockNames : str
            Names of each mesh block
        cellBounds : 
        elementTypes : 
            Element types of each mesh block
        numberOfCells : int
            Total number of cells
        numberOfPoints : int
            Total number of points
        fieldSpecifications : dict
            Dict containing the mesh field specifications
    """
    def __init__(self, xml):
        """
        Parameters
        ----------
            xml : XML
                XML object containing the information on the mesh
        """
        self.xml = xml
        
        mesh = xml.mesh["InternalMesh"]
        elementRegion = xml.elementRegions["CellElementRegion"]
        fieldSpecifications = xml.fieldSpecifications

        name = mesh["name"]

        xCoords = list(eval(mesh["xCoords"]))
        yCoords = list(eval(mesh["yCoords"]))
        zCoords = list(eval(mesh["zCoords"]))

        self.bounds = [[xCoords[0], xCoords[-1]], [yCoords[0], yCoords[-1]], [zCoords[0], zCoords[-1]]]


        nxStr = mesh["nx"].strip('').replace('{','').replace('}','').split(',')
        nyStr = mesh["ny"].strip('').replace('{','').replace('}','').split(',')
        nzStr = mesh["nz"].strip('').replace('{','').replace('}','').split(',')

        nx = [eval(nx) for nx in nxStr]
        ny = [eval(ny) for ny in nyStr]
        nz = [eval(nz) for nz in nzStr]

        self.nx = nx
        self.ny = ny
        self.nz = nz

        order = 1
        self.order = order

        self.cellBlockNames = mesh["cellBlockNames"].strip('').replace('{','').replace('}','').split(',')

        xlayers = []
        ylayers = []
        zlayers = []
        for i in range(len(nx)):
            xlayers.append([xCoords[i], xCoords[i+1]])
        for i in range(len(ny)):
            ylayers.append([yCoords[i], yCoords[i+1]])
        for i in range(len(nz)):
            zlayers.append([zCoords[i], zCoords[i+1]])

        self.layers = [xlayers, ylayers, zlayers]

        xCellsBounds = np.zeros(sum(nx)+1)
        yCellsBounds = np.zeros(sum(ny)+1)
        zCellsBounds = np.zeros(sum(nz)+1)

        for i in range(len(nx)):
            xstep = (xlayers[i][1]-xlayers[i][0])/nx[i]
            if i == 0:
                xCellsBounds[0:nx[i]] = np.arange(xlayers[i][0], xlayers[i][1], xstep)
            else :
                xCellsBounds[nx[i-1]:sum(nx[0:i+1])] = np.arange(xlayers[i][0], xlayers[i][1], xstep)
        xCellsBounds[nx[-1]] = xlayers[i][1]

        for i in range(len(ny)):
            ystep = (ylayers[i][1]-ylayers[i][0])/ny[i]
            if i == 0:
                yCellsBounds[0:ny[i]] = np.arange(ylayers[i][0], ylayers[i][1], ystep)
            else :
                xCellsBounds[ny[i-1]:sum(ny[0:i+1])] = np.arange(ylayers[i][0], ylayers[i][1], ystep)
        yCellsBounds[ny[-1]] = ylayers[i][1]

        for i in range(len(nz)):
            zstep = (zlayers[i][1]-zlayers[i][0])/nz[i]
            if i == 0:
                zCellsBounds[0:nz[i]] = np.arange(zlayers[i][0], zlayers[i][1], zstep)
            else :
                zCellsBounds[nz[i-1]:sum(nz[0:i+1])] = np.arange(zlayers[i][0], zlayers[i][1], zstep)
        zCellsBounds[nz[-1]] = zlayers[i][1]

        self.cellBounds = [xCellsBounds, yCellsBounds, zCellsBounds]

        elementTypes = mesh["elementTypes"].strip('').replace('{','').replace('}','').split(',')

        self.elementTypes=[]
        for type in elementTypes:
            if type == "C3D8":
                self.elementTypes.append("Hexahedron")
            else:
                self.elementTypes.append(type)


        self.numberOfCells = sum(nx) * sum(ny) * sum(nz)
        self.numberOfPoints = (sum(nx) + 1) * (sum(ny) + 1) * (sum(nz) + 1)

        self.fieldSpecifications = xml.fieldSpecifications
        self.isSet = True


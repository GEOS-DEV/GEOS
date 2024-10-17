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

import os
import vtk
from vtk.util import numpy_support as VN

from .VtkFieldSpecifications import VTKCellSpecifications, VTKPointSpecifications
from utilities.model.utils.vtkUtils import cGlobalIds



class VTKMesh:
    """
    VTK format Mesh. Now handling .vtu .vts .pvts .pvtu

    Attributes
    ----------
        meshfile : str
            Mesh filename
        vtktype : str
            Format of the VTK mesh
        bounds : tuple of int 
            Real bounds of the mesh (xmin, xmax, ymin, ymax, zmin, zmax)
        numberOfPoints : int
            Total number of points of the mesh
        numberOfCells : int
            Total number of cells of the mesh
        isSet : bool
            Whether or not the mesh properties have been set
        hasLocator : bool
            Whether or not the mesh cell locator has been initialized
    """
    def __init__(self, meshfile):
        """
        Parameters
        ----------
            meshfile : str
                Mesh filename
        """    
        self.meshfile = meshfile
        self.vtktype = os.path.splitext(self.meshfile)[-1][1:]

        self.bounds = None
        self.numberOfPoints = None
        self.numberOfCells = None

        self.isSet = False
        self.hasLocator = False


    def getReader(self):
        """Return the appropriate reader given the VTK format
        
        Returns
        --------
            vtk.vtkXMLReader
                Appropriate reader given the format of the VTK mesh file.
        """
        
        if self.vtktype == "vtu":
                return vtk.vtkXMLUnstructuredGridReader()
        elif self.vtktype == "vts":
                return vtk.vtkXMLStructuredGridReader()
        elif self.vtktype == "pvtu":
                return vtk.vtkXMLPUnstructuredGridReader()
        elif self.vtktype == "pvts":
                return vtk.vtkXMLPStructuredGridReader()  
        else:
            print("This VTK format is not handled.")
            return None

    
    def read(self):
        """Read information from the VTK file

        Returns
        --------
            vtk.vtkDataObject
                General representation of VTK mesh data 
        """
        reader = self.getReader()
        reader.SetFileName(self.meshfile)
        reader.Update()

        return reader.GetOutput()


    def setMeshProperties(self):
        """Read and set as attributes the bounds, number of points and cells"""
        data = self.read()

        self.bounds = data.GetBounds()
        self.numberOfPoints = data.GetNumberOfPoints()
        self.numberOfCells = data.GetNumberOfCells()

        self.isSet = True
  
 
    def getBounds(self):
        """
        Get the bounds of the mesh in the format:
            (xmin, xmax, ymin, ymax, zmin, zmax)

        Returns
        ------- 
            tuple or None
                Bounds of the mesh
        """
        return self.bounds


    def getNumberOfPoints(self):
        """
        Get the total number of points of the mesh

        Returns
        -------
            int
                Number of points
        """
        return self.numberOfPoints


    def getNumberOfCells(self):
        """
        Get the total number of cells of the mesh

        Returns
        -------
            int
                Number of cells
        """
        return self.numberOfCells 


    def getCellData(self):
        """Read the cell data
    
        Returns
        --------
            VTKCellSpecifications
                Cell data information
        """
        data = self.read()

        return VTKCellSpecifications(data.GetCellData())


    def getPointData(self):
        """Read the point data

        Returns
        --------
            VTKPointSpecifications
                Point data information
        """
        data = self.read()

        return VTKPointSpecifications(data.GetPointData())
    

    def getArray(self, name, dtype="cell", copy=False, sorted=False):
        """
        Return a cell or point data array. If the file is a pvtu, the array is sorted with global ids

        Parameters
        -----------
            name : str
                Name of the vtk cell/point data array
            dtype : str
                Type of vtk data \
                `cell` or `point`
            copy : bool
                Return a copy of the requested array
                Default is False
            sorted : bool
                Return the array sorted with respect to GlobalPointIds or GlobalCellIds
                Default is False

        Returns
        --------
            array : numpy array
                Requested array
        """
        assert dtype.lower() in ("cell", "point")

        if dtype.lower() == "cell":
            fdata = self.getCellData()
        else:
            fdata = self.getPointData()

        if copy:
            array = fdata.getCopyArray(name, sorted=sorted)
        else:
            array = fdata.getArray(name, sorted=sorted)

        return array


    def extractMesh(self, center, srootname, dist=[None, None, None], comm=None, export=True):
        """
        Extract a rectangular submesh such that for each axis we have the subax: [center-dist, center+dist]

        Parameters
        ---------
            center : 3d float
                Requested center of the subbox 
            srootname : str
                Submesh root filename
            dist : 3d float
                Distance to the center in each direction
            comm : MPI.COMM_WORLD
                MPI communicator

        Returns
        -------
            VTKMesh 
                Submesh extracted
        """
        assert self.vtktype in ("vtu", "pvtu", "vts", "pvts")
        vtktype = self.vtktype[-3:]
        sfilename = ".".join((srootname, vtktype))

        if comm is None or comm.Get_rank() == 0:
            if not self.isSet:
                self.setMeshProperties()
            minpos = []
            maxpos = []
            
            for i in range(3):
                xmin, d = self.getSubAx(center[i], dist[i], ax=i+1)
                minpos.append(xmin)
                maxpos.append(xmin+d)

            data = self.read()
            submesh = VTKSubMesh(sfilename, data, minpos, maxpos, create=export)

        else:
            submesh = VTKMesh(sfilename)

        # if file creation, wait for rank 0 to finish
        if export:
            info = "Done"
            comm.bcast(info, root=0)

        return submesh
   
 
    def getSubAx(self, center, dist, ax):
        """
        Return the min and max positions in the mesh given the center, distance and ax considered. If the 2*distance if greater than the bounds, the min/max is the corresponding mesh bound.

        Parameters
        ----------
            center : float
                Central position considered
            dist : float
                Max distance requested
            ax : int
                Ax to consider (1, 2, 3)
        
        Returns
        -------
            min, max : float
                Min and Max positions 
        """
        assert(type(ax) == int)
   
        bounds = [self.bounds[(ax-1)*2], self.bounds[ax*2-1]]

        if dist is not None:
            dist = abs(dist)
            ox = max(bounds[0], center-dist)
            x = min(bounds[1]-ox, 2*dist)
        else:
            ox = bounds[0]
            x = bounds[1]

        return ox, x


    def getNumberOfBlocks(self):
        """Return the number of blocks of a mesh."""
        if self.vtktype in ["pvtu", "pvts"]:
            with open(self.meshfile) as ff: 
                nb = 0 
                for line in ff: 
                    m = line.split()
                    if m[0] == '<Piece':
                        nb+=1
            return nb

        else:
            return 1


    def getGlobalIds(self, dtype="cell"):
        """Return the global ids of the cells or points. If the mesh is an extract of an original mesh, it is the local to global map

        Parameters
        ----------
            dtype : str
                Type of data: `cell` or `point`

        Returns
        --------
            array-like
                Global Ids
        """
        assert dtype.lower() in ("cell", "point")

        if dtype.lower() == "cell":
            fdata = self.getCellData()
        else:
            fdata = self.getPointData()
            
        if fdata.hasArray(f"Global{dtype.title()}Ids"):
            gids = fdata.getArray(f"Global{dtype.title()}Ids").ravel()
            return gids
        
        else:
            print("No global Ids array found in this VTK mesh")
            return None


    def getExtractToGlobalMap(self):
        """Return the global cell ids
        
        Returns
        --------
            array-like : 
                Global cell Ids or None if not set in the mesh
        """
        return self.getGlobalIds()


    def export(self, data=None, rootname=None, vtktype=None):
        """
        Write VTK data in a file

        Parameters
        ----------
            data : vtkDataSet
                vtk.vtkStructuredGrid or vtk.vtkUnstructuredGrid
                Default is self.read()
            rootname : str
                Root of the output filename
                Default is self.meshfile (without extension)
            vtktype : str
                Format of the output VTK
                Default is self.vtktype

        Returns
        --------
            filename : str
                Output filename
        """
        if vtktype is None:
            vtktype = self.vtktype
        if rootname is None:
            rootname, _ = os.path.splitext(self.meshfile)
        if data is None:
            data = self.read()
        
        filename = ".".join((rootname, vtktype))

        writer = self.getWriter(vtktype=vtktype)
        writer.SetFileName(filename)
        writer.SetInputData(data)
        writer.Update()
        writer.Write()

        return filename


    def getWriter(self, vtktype=None):
        """
        Return the VTK writer

        Returns
        --------
            vtk.vtkXMLWriter
                Appropriate writer given the format of the VTK mesh file.
        """
        if vtktype is None:
            vtktype = self.vtktype

        if vtktype == "vts":
            return vtk.vtkXMLStructuredGridWriter()
        elif vtktype == "vtu":
            return vtk.vtkXMLUnstructuredGridWriter()


    def setCellLocator(self):
        """Set the cell locator"""
        if not self.isSet:
          self.setMeshProperties()

        if not self.hasLocator:
          self.cellLocator = vtk.vtkCellLocator()
          self.cellLocator.SetDataSet(self.read())
          self.cellLocator.BuildLocator()
          self.hasLocator = True


    def getCellContainingPoint(self, point):
        """
        Return the global index of the cell containing the coordinates

        Parameters
        -----------
            point : array-like of float
                Point coordinates

        Returns
        --------
            cellIds : int
                id of the cell containing the given point
        """
        if not self.hasLocator:
          self.setCellLocator()

        cellIds = self.cellLocator.FindCell([point[0], point[1], point[2]])
        return cellIds


    def interpolateValues(self, centers, name, values):
        """
        Interpolate the given cell data over the given points

        Parameters
        -----------
            centers : list of list of float
                Center coordinates
            name : str
                Name of the new array
            values : numpy array
                New values
        
        Returns
        --------
           interpValues : 
                interpolated values over the given points 
        """
        if not self.isSet:
            self.setMeshProperties()

        dest = vtk.vtkPointSet()
        destPoints = vtk.vtkPoints()

        for point in centers:
            destPoints.InsertNextPoint([point[0], point[1], point[2]])
        dest.SetPoints(destPoints)

        transferArray = vtk.vtkDoubleArray()
        for value in values:
            transferArray.InsertNextTuple1(value)
        transferArray.SetName(name)

        data = self.read()
        data.GetCellData().AddArray(transferArray)
        resample = vtk.vtkResampleWithDataSet()
        resample.SetSourceData(data)
        resample.SetInputData(dest)
        resample.Update()
        
        pointdata = resample.GetOutput().GetPointData()

        interpValues = None
        for i in range(pointdata.GetNumberOfArrays()):
            array = pointdata.GetArray(i)
            if array.GetName() == name:
                interpValues = VN.vtk_to_numpy(array)

        return interpValues



class VTKSubMesh(VTKMesh):
    """
    Class defining a submesh of an existing VTK mesh

    Attributes
    -----------
        meshfile : str
            Submesh filename
        vtktype : str
            Format of the VTK submesh
        bounds : tuple of int 
            Real bounds of the mesh (xmin, xmax, ymin, ymax, zmin, zmax)
        numberOfPoints : int
            Total number of points of the submesh
        numberOfCells : int
            Total number of cells of the submesh
        isSet : bool
            Whether or not the mesh properties have been set
    """
    def __init__(self, meshfile, data, minpos, maxpos, create=True):
        """
        Parameters
        -----------
            meshfile : str
                Submesh filename
            data : vtk.vtkDataObject
                General representation of the original mesh
            minpos : 3d array-like of float
                Minimal positions for the cropping for each axis
            maxpos : 3d array-like of float
                Maximal positions for the cropping for each axis
            create : bool
                Whether or not to create the VTKfile
                Default is True
        """
        super().__init__(meshfile)

        sdata = self.__setData(data, minpos, maxpos)
        self.__setGlobalIds(sdata, data)

        if create:
            self.export(data=sdata)


    def __setGlobalIds(self, sdata, data):
        """
        Set the global cell Ids of the submesh

        Parameters
            sdata : vtk.vtkDataObject
                General representation of the submesh
        """
        if self.vtktype == "vtu":
            subcdata = sdata.GetCellData()
            if subcdata.HasArray("GlobalCellIds") == 1:
                subcdata.RemoveArray("vtkOriginalCellIds")
            else:
                cgids = subcdata.GetArray("vtkOriginalCellIds")
                cgids.SetName("GlobalCellIds")

        elif self.vtktype == "vts":
            if not sdata.GetCellData().HasArray("GlobalCellIds"):
                nx_extract, ny_extract, nz_extract = sdata.GetDimensions()
                dx = data.GetBounds()[1]/data.GetExtent()[1]
                dy = data.GetBounds()[3]/data.GetExtent()[3]
                dz = data.GetBounds()[5]/data.GetExtent()[5] 
                nx, ny, nz = data.GetDimensions()
                xmin, ymin, zmin = sdata.GetBounds()[0::2]
                xmin0, ymin0, zmin0 = self.bounds[0::2]
                
                cgids = cGlobalIds(nx_extract, ny_extract, nz_extract,
                                    dx, dy, dz,
                                    xmin, ymin, zmin,
                                    nx, ny, nz,
                                    xmin0, ymin0, zmin0)

                subcdata = sdata.GetCellData()
                cgidsAsVtkArray = VN.numpy_to_vtk(num_array=cgids.ravel(), deep=True)
                cgidsAsVtkArray.SetName("GlobalCellIds")
                subcdata.AddArray(cgidsAsVtkArray)
    

    def __setData(self, data, minpos, maxpos):
        """
        Return the submesh extracted from the whole mesh dataset

        Parameters
        -----------
            data : vtk.vtkDataObject
                General representation of the original mesh
            minpos : 3d array-like of float
                Minimal positions for the cropping for each axis
            maxpos : 3d array-like of float
                Maximal positions for the cropping for each axis
        """
        assert None not in minpos and len(minpos) == 3
        assert None not in maxpos and len(maxpos) == 3

        if self.vtktype == "vtu":
            cellLocator = vtk.vtkCellLocator()
            cellLocator.SetDataSet(data)
            cellLocator.BuildLocator()

            idList = vtk.vtkIdList()
            cells = cellLocator.FindCellsWithinBounds([minpos[0], maxpos[0], minpos[1], maxpos[1], minpos[2], maxpos[2]], idList)

            
            #Extraction of the cells
            extract = vtk.vtkExtractCells()
            extract.SetInputData(data)
            extract.SetCellList(idList)
            extract.Update()

            dataExtract = extract.GetOutput()
            

        elif self.vtktype == "vts":
            # vtkExtractGrid requires the [i,j,k] coordinates
            # distances and positions have to be converted 
            dx = data.GetBounds()[1]/data.GetExtent()[1]
            dy = data.GetBounds()[3]/data.GetExtent()[3]
            dz = data.GetBounds()[5]/data.GetExtent()[5] 

            minx = int(minpos[0]//dx)
            miny = int(minpos[1]//dy)
            minz = int(minpos[2]//dz)

            maxx = int(maxpos[0]//dx)
            maxy = int(maxpos[1]//dy)
            maxz = int(maxpos[2]//dz)

            #Extraction of the grid
            extract = vtk.vtkExtractGrid()
            extract.SetInputData(data)

            extract.SetVOI(minx, maxx, miny, maxy, minz, maxz)
            extract.Update()

            dataExtract = extract.GetOutput()

        return dataExtract
    

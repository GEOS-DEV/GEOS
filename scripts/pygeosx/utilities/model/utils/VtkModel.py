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
import sys
import numpy as np
import vtk
from . import vtkUtils as vtkwriter


class VTKModel:
    """
    Define a VTK model

    Attributes
    -----------
        filename : str
            Name of the VTK file with extension
        rootname : str
            Root of the filename
        vtktype : str
            Type of VTK format
        cellData : dict
            Contains all cell data arrays
        pointData : dict
            Contains all point data arrays     
        cgids : list of int
            Cell global Ids
        pgids : list of int
            Point global Ids
        bmin : tuple of float
            Bound min for each dimension
        glmin : tuple of float
            Bound min for each dimension for the global model
        n : tuple of int
            Number of elements for each dimension
        gln : tuple of int
            Number of elements for each dimension for the global model
        d : tuple of float
            Step size for each dimension
        vertices : tuple of arrays
            Model points coordinates
        ctype : array-like
            Cells datatype
    """
    def __init__(self, vtkfile, cellData=None, pointData=None):
        """
        Parameters
        ----------
            vtkfile : str
                Filename
            cellData : dict
                Dictionary containing the cell data
            pointData : dict
                Dictionary containing the point data
        """
        self.filename = vtkfile
        self.rootname, ext = os.path.splitext(self.filename) 
        self.vtktype = ext[1:]
        self.cellData = cellData
        self.pointData = pointData
        self.cgids = None
        self.pgids = None

        self.bmin = (None, None, None)
        self.n = (None, None, None)
        self.d = (None, None, None)
        self.gln = (None, None, None)
        self.glmin = (None, None, None)   
        self.vertices = None


    def getReader(self):
        """
        Returns the vtkXMLReader corresponding to the VTK format
        
        Returns
        -------
            vtkXMLReader
                Reader for this VTK format
        """
        if hasattr(self, "reader"):
            return self.reader
        else:
            if self.vtktype == "vtu":
                    return vtk.vtkXMLUnstructuredGridReader()
            elif self.vtktype == "vts":
                    return vtk.vtkXMLStructuredGridReader()
            else:
                print("This VTK format is not handled.")
                return None


    def setNumberOfElements(self, n):
        """
        Define the number of points, cells, and total number of cells

        Parameters
        ----------
            n : array-like of int
                Number of points
        """
        self.n = n
        self.nc = tuple(np.array(n)-1)
        self.ncells = np.prod(self.nc)


    def setGlobalNumberOfElements(self, gln):
        """
        Define the number of elements of the global model

        Parameters
        ----------
            n : array-like of int
                Number of points
        """
        self.gln = gln
    

    def setOrigin(self, bmin):
        """
        Define the origin of the model

        Parameters
        -----------
            bmin : array-like of float
                Origin of the model
        """
        self.bmin = bmin


    def setStepSize(self, d):
        """
        Define the step size of the model

        Parameters
        -----------
            d : array-like of float
                Step sizes
        """
        self.d = d


    def setGlobalOrigin(self, glmin):
        """
        Define the origin of the global model

        Parameters
        -----------
            glmin : array-like of float
                Global origin of the model
        """
        self.glmin = glmin


    def getIndexMin(self):
        """
        Return the minimal indices of the model

        Returns
        --------
            array-like of int
        """
        if None not in self.d and None not in self.bmin:
            bmin = np.array(self.bmin)
            d = np.array(self.d)

            imin = np.array(bmin//d, dtype=int)
            return tuple(imin)
        else:
            return (None, None, None)


    def getIndexMax(self):
        """
        Return the maximal indices of the model

        Returns
        --------
            array-like of int
        """
        if None not in self.d and None not in self.bmin and None not in self.n:
            bmin = np.array(self.bmin)
            d = np.array(self.d)
            n = np.array(self.n)

            imax = np.array(bmin//d, dtype=int)+(n-1)
            return tuple(imax)
        else:
            return (None, None, None)


    def setVertices(self):
        """
        Compute point coordinates (vertices)

        Parameters
        ----------
            n : array-like of int
                Number of elements in each dimension
            d : array-like of float
                Step sizes for each dimension
            bmin : array-like of float
                Minimal bound for each dimension
        """
        nx, ny, nz = self.n  
        dx, dy, dz = self.d
        xmin, ymin, zmin = self.bmin

        if self.vtktype == "vts":
            v1, v2, v3 = vtkwriter.x_y_z(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin)
        elif self.vtktype == 'vtu':
            v1, v2, v3 = vtkwriter.xyz(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin)

        self.vertices = (v1, v2, v3)
    

    def setCellGlobalIds(self):
        """
        Compute the global ids of the cells

        Parameters
        ----------
            n : array-like of int
                Number of elements in each dimension
            d : array-like of float
                Step sizes for each dimension
            bmin : array-like of float
                Minimal bound for each dimension
            gln : array-like of int
                Global number of elements for each dimension
            glmin : array-like of float
                Global minimal bound for each dimension
        """
        nx, ny, nz = self.n  
        dx, dy, dz = self.d
        xmin, ymin, zmin = self.bmin
        glNx, glNy, glNz = self.gln
        xmin0, ymin0, zmin0 = self.glmin

        assert None not in self.n
        assert None not in self.d
        assert None not in self.bmin
        assert None not in self.gln
        assert None not in self.glmin


        self.cgids = vtkwriter.cGlobalIds(nx, ny, nz, dx,dy,dz,xmin,ymin,zmin,glNx,glNy, glNz,xmin0,ymin0,zmin0)

        if not self.cellData:
            self.cellData = {"GlobalCellIds" : self.cgids}
        else:
            self.cellData = {**self.cellData, **{"GlobalCellIds" : self.cgids}}


    def getCellGlobalIds(self):
        """
        Return the value of the cell global Ids

        Returns
        -------
            array-like of int
        """
        return self.cgids


    def setPointGlobalIds(self):
        """
        Compute the global ids of the points

        Parameters
        ----------
            n : array-like of int
                Number of elements in each dimension
            d : array-like of float
                Step sizes for each dimension
            bmin : array-like of float
                Minimal bound for each dimension
            gln : array-like of int
                Global number of elements for each dimension
            glmin : array-like of float
                Global minimal bound for each dimension
        """
        nx, ny, nz = self.n  
        dx, dy, dz = self.d
        xmin, ymin, zmin = self.bmin
        glNx, glNy, glNz = self.gln
        xmin0, ymin0, zmin0 = self.glmin

        self.pgids = vtkwriter.pGlobalIds(nx, ny, nz, dx,dy,dz,xmin,ymin,zmin,glNx,glNy, glNz,xmin0,ymin0,zmin0)

        if not self.pointData:
            self.pointData = {"GlobalPointIds" : self.pgids}
        else:
            self.pointData = {**self.pointData, **{"GlobalPointIds" : self.pgids}}
    

    def getPointGlobalIds(self):
        """
        Return the value of the points global Ids

        Returns
        -------
            array-like of int
        """
        return self.pgids


    def addCellData(self, ckey, cellData):
        """
        Append or update cell data
 
        Parameters
        ----------
            ckey : str
                Cell data name
            cellData : array-like
                Cell data array
        """
        if not self.cellData:
            self.cellData = {ckey : cellData}
        else:
            self.cellData = {**self.cellData, ckey : cellData}


    def addPointData(self, pkey, pointData):
        """
        Append or update point data
 
        Parameters
        ----------
            pkey : str
                Point data name
            pointData : array-like
                Point data array
        """
        if not self.pointData:
            self.pointData = {pkey : pointData}
        else:
            self.pointData = {**self.pointData, pkey : pointData}


    def setCellsType(self, ncells, ctype=12):
        """
        Set cell type of all the cells

        Parameters  
        -----------
            ncells : int
                Number of cells
            ctype : int
                Type of the cells 
                12 for vtk.VtkHexahedron.tid
        """
        self.ctype = np.full((ncells), fill_value=ctype, dtype='uint8')


    def export(self):
        pass        



class VTUModel(VTKModel):
    """
    Define a structured grid model
    Inheritance from VTKModel

    Attributes
    -----------
        filename : str
            Name of the VTK file
        rootname : str
            Root of the filename
        vtktype : str
            Type of VTK format (unstructured grid = "vtu")
        cellData : dict
            Contains all cell data arrays
        pointData : dict
            Contains all point data arrays 
        reader : vtkXMLUnStructuredGridReader
            VTK file reader     
        cgids : list of int
            Cell global Ids
        pgids : list of int
            Point global Ids
        bmin : tuple of float
            Bound min for each dimension
        glmin : tuple of float
            Bound min for each dimension for the global model
        n : tuple of int
            Number of elements for each dimension
        gln : tuple of int
            Number of elements for each dimension for the global model
        d : tuple of float
            Step size for each dimension
        vertices : tuple of arrays
            Model points coordinates
        ctype : array-like
            Cells datatype
        connectivity : array-like
            Connectivities
        offsets : array-like
            Cells offsets 
    """
    def __init__(self, vtkfile, cellData=None, pointData=None):
        """
        Parameters
        ----------
            vtkfile : str
                Filename
            cellData : dict (optional)
                Contains all cell data arrays
            pointData : dict (optional)
                Contains all point data arrays
        """
        VTKModel.__init__(self, vtkfile, cellData, pointData)
        self.vtktype = "vtu"
        self.reader = vtk.vtkXMLUnstructuredGridReader()

        self.connectivity = None
        self.offsets = None


    def setConnectivity(self):
        """
        Compute the connectivity and offsets

        Parameters
        -----------
            n : array-like of int
                Number of elements for each dimension

        Returns
        --------
            conn : array-like 
                Connectivities
            offset : array-like
                Offsets
            numCells : int
                Number of cells
        """
        nx, ny, nz = self.n
        conn, offset, ncells = vtkwriter.connectivity(nx, ny,nz)
        
        self.connectivity = conn
        self.offsets = offset
        self.ncells = ncells


    def setVertices(self):
        """
        Compute point coordinates (vertices) for a vtu format

        Parameters
        ----------
            n : array-like of int
                Number of elements in each dimension
            d : array-like of float
                Step sizes for each dimension
            bmin : array-like of float
                Minimal bound for each dimension
        """
        nx, ny, nz = self.n
        dx,dy,dz = self.d
        xmin,ymin,zmin = self.bmin

        x,y,z = vtkwriter.xyz(nx, ny, nz, dx, dy, dz, xmin, ymin, zmin)

        self.vertices = (x,y,z)
        
    
    def export(self, gids=False):
        """
        Write the .vtu file with celldata and pointdata defined in the class

        Parameters
        ----------
            vertices : array-like
                Coordinates of the points along each dimension
            connectivity : array-like
                Cells point connectivity
            offsets : array-like
                Offset into the connectivity of the cells
            pgids : array-like
                Point global Ids
            cgids : array-like
                Cell global Ids
        """
        self.setVertices()
        x,y,z = self.vertices

        self.setConnectivity()
                
        if gids:
            self.setCellGlobalIds()
            self.setPointGlobalIds()

        vtkwriter.unstructuredGridToVTK(self.rootname, x, y, z, connectivity=self.connectivity, offsets=self.offsets, cell_types=self.ctype, pointData=self.pointData, cellData=self.cellData)



class VTSModel(VTKModel):
    """
    Define a structured grid model
    Inheritance from VTKModel

    Attributes
    -----------
        filename : str
            Name of the VTK file
        rootname : str
            Root of the filename
        vtktype : str
            Type of VTK format (structured grid = "vts")
        cellData : dict
            Contains all cell data arrays
        pointData : dict
            Contains all point data arrays 
        reader : vtkXMLUnStructuredGridReader
            VTK file reader  
        cgids : list of int
            Cell global Ids
        pgids : list of int
            Point global Ids
        bmin : tuple of float
            Bound min for each dimension
        glmin : tuple of float
            Bound min for each dimension for the global model
        n : tuple of int
            Number of elements for each dimension
        gln : tuple of int
            Number of elements for each dimension for the global model
        d : tuple of float
            Step size for each dimension
        vertices : tuple of arrays
            Model points coordinates    
    """
    def __init__(self, vtkfile, cellData=None, pointData=None):
        """
        Parameters
        ----------
            vtkfile : str
                Filename
            cellData : dict
                Contains all cell data arrays
            pointData : dict
                Contains all point data arrays
        """
        VTKModel.__init__(self, vtkfile, cellData, pointData)
        self.vtktype = "vts"
        self.reader = vtk.vtkXMLStructuredGridReader()


    def setVertices(self):
        """
        Compute point coordinates (vertices) for a vts format

        Parameters
        ----------
            n : array-like of int
                Number of elements in each dimension
            d : array-like of float
                Step sizes for each dimension
            bmin : array-like of float
                Minimal bound for each dimension

        Returns
        -------
            3d array-like
                Arrays of point coordinates along each dimension
        """
        nx, ny, nz = self.n
        dx,dy,dz = self.d
        xmin,ymin,zmin = self.bmin

        x,y,z = vtkwriter.x_y_z(nx,ny,nz, dx,dy,dz,xmin,ymin,zmin)
        self.vertices = (x,y,z)


    def export(self, gids=False):
        """
        Write the .vts file with celldata and pointdata defined in the class

        Parameters
        ----------
            vertices : array-like
                Coordinates of the points along each dimension
            pgids : array-like
                Point global Ids
            cgids : array-like
                Cell global Ids
        """
        self.setVertices()
        x,y,z = self.vertices
        
        if gids:
            self.setCellGlobalIds()
            self.setPointGlobalIds()
        
        imin = self.getIndexMin()
        if None in imin:
            print("Problem in the determination of the block index origin.")
            sys.exit(1)

        vtkwriter.structuredToVTK(self.rootname, x, y, z, cellData=self.cellData, pointData=self.pointData, start=imin)



class PVTKModel:
    """
    Parallel VTK model

    Attributes
    ----------
        filename : str
            Filename of the pvtk header
        rootname : str
            Root of the filename of all files from the PVTK model
        pvtktype : str
            VTK type ('pvtu' or 'pvts')
        vtkfiles : dict
            Contains all sources filenames
            Block number ID out of the total number of blocks (key) associated to the corresponding VTK filename
        vtkmodels : dict
            Contains all sources VTKModels
            Block number ID out of the total number of blocks (key) associated to corresponding VTKModel (VTUModel or VTSModel)
        cellInfo : dict
            Contains cell data information
            Cell array names (key) associated to data type and number of components
        pointInfo : dict
            Contains point data information
            Point array names (key) associated to data type and number of components
        starts : dict
            Contains the blocks min indices
            Block number IDs for each dimension (key) associated to min indices
        ends : dict
            Contains the blocks max indices
            Block number IDs for each dimension (key) associated to max indices
        nblocks : tuple of int
            Number of blocks for each dimension
        nbltot : int
            Total number of blocks
        n : tuple of int
            Whole extent of the model
    """
    def __init__(self, pvtkfile, vtkfiles=None):
        """
        Parameters
        ----------
            pvtkfile : str
                Filename of the PVTK model
            vtkfiles : array-like of str (optional)
                List of filenames of VTK files constituting the PVTK model
        """
        self.filename = pvtkfile
        self.rootname, ext = os.path.splitext(self.filename)
        self.pvtktype = ext[1:]

        self.vtkfiles = None
        self.__setSources(vtkfiles)

        self.vtkmodels = None
        self.cellInfo = None
        self.pointInfo = None
        self.starts = None
        self.ends = None
        


    def __setSources(self, vtkfiles=None):
        """
        Set the VTK files associated to this PVTK

        Parameters
        ----------
            vtkfiles : array-like of str
                Sorted list of files that constitute the PVTK
        """
        if not vtkfiles:
            try:
                vtkfiles = self.read()
            except:
                pass
        
        if vtkfiles:
            self.vtkfiles = {id: vtkf for id, vtkf in enumerate(vtkfiles)}


    def __setVtkModels(self):
        """
        Set the VTK models associated to this PVTK from the vtkfiles attribute

        Only useful if we need to export all the blocks from the PVTKModel
        """
        if self.vtkfiles:
            if self.pvtktype == "pvtu":
                vtkmodels = {id: VTUModel(vtkfile) for id, vtkfile in self.vtkfiles.items()}
            else:
                vtkmodels = {id: VTSModel(vtkfile) for id, vtkfile in self.vtkfiles.items()}

        self.vtkmodels = {**self.vtkmodels, **vtkmodels}
        
        if self.vtkmodels:
            for id, model in self.vtkmodels.items():
                self.addBlockInfo(block=model, bn=id)
        

    def setNumberOfBlocks(self, nblocks):
        """
        Set the number of blocks for each dimension and the total number of blocks of the model

        Parameters
        ----------
            nblocks : array-like of int
                Number of blocks for each dimension
        """
        self.nblocks = nblocks
        self.nbltot = np.prod(nblocks)
        

    def setWholeExtent(self, extent):
        """
        Define the extent of the global model

        Parameters
        -----------
            extent : array-like of int
                Min and Max indices of the global model
        """
        self.n = extent


    def addCellInfo(self, cinfo=None):
        """
        Add cell data name and type to the cell information attribute

        Parameters
        ----------
            cinfo : dict, optional
                Cell array name as key, cell datatype and number of component as value
        """
        if cinfo:
            if self.cellInfo:
                for k, inf in cinfo.items():
                    if k in self.cellInfo.keys():
                        datatype, ncomp = inf 
                        dref, ncref = self.cellInfo[k] 
                        assert datatype == dref
                        assert ncomp == ncref
                    else:
                        self.cellInfo = {**self.cellInfo, k: inf}
            else:
                self.cellInfo = cinfo


    def addPointInfo(self, pinfo=None):
        """
        Add point data name and type to the point information attribute

        Parameters
        ----------
            pinfo : dict, optional
                Point array name as key, data type and number of components as value
        """
        if pinfo:
            if self.pointInfo:
                for k, inf in pinfo.items():
                    if k in self.pointInfo.keys():
                        datatype, ncomp = inf 
                        dref, ncref = self.pointInfo[k] 
                        assert datatype == dref
                        assert ncomp == ncref
                    else:
                        self.cellInfo = {**self.pointInfo, k: inf}
            else:
                self.pointInfo = pinfo

    
    def _addBlockExtents(self, bijk, start, end):
        """
        Add the block minimal and maximal indices of a block source

        Parameters
        -----------
            bijk : tuple of int 
                Block number IDs for each dimension
            start : tuple of int
                Block minimal indices
            end : tuple of int
                Block maximal indices
        """
        if self.pvtktype == "pvts":
            if not self.starts or not self.ends:
                self.starts = {}
                self.ends = {}

            self.starts = {**self.starts, (*bijk, 0): start}
            self.ends = {**self.ends, (*bijk, 1): end}


    def readHeader(self):
        """
        Return the PVTK file header

        Returns
        -------
            str
        """
        with open(self.filename, "r") as f:
            ffull = f.read()
        return ffull


    def getSourceFilesFromHeader(self):
        """
        Return the source files from the PVTK file header

        Returns
        -------
            list of str
        """
        pvtk = self.readHeader()
        vtkfiles = []

        for l in pvtk.split("\n"):
            if "Source=" in l:
                _, vtkfile, _ = l.split('\"')
                vtkfiles += [vtkfile]   

        return vtkfiles


    def addBlockModel(self, block, bn):
        """
        Add a VTK model source to the PVTK model

        Parameters
        ----------
            block : VTKModel
                VTKModel of a specific block from the global PVTK model
            bn : int
                Block number Id out of the total number of blocks
        """
        if not self.vtkmodels:
            self.vtkmodels = {bn: block}
        else:
            self.vtkmodels = {**self.vtkmodels, **{bn: block}}

        self.addBlockInfo(block, bn=bn)


    def addBlockInfo(self, block, bijk=None, bn=None):
        """
        Add information about a specific block (cell and point informations, min and max indices)

        Parameters
        -----------
            block : VTKModel
                VTKModel of a specific block from the global PVTK model
            bijk : tuple of int (optional)
                Block number Ids for each dimension
                Required if bn is not given
            bn : int (optional)
                Block number Id out of the total number of blocks
                Required if bijk is not given
        """
        if not bijk and not bn:
            raise ValueError("The block identification is required. You can set bijk = (ni, nj, nk) the block ids for each dimension, or bn = int, the id on the total number of blocks")

        if not bijk and bn:
            nijk = self.getAllBlocksIndices()
            bijk = nijk[bn]
        
        else:
            nijk = self.getAllBlocksIndices()
            bn = nijk.index(bijk)

        cinfo = {}
        pinfo = {}

        if block.cellData:
            for kc, cdata in block.cellData.items():
                    cinfo = {**cinfo, kc : (cdata.dtype, len(np.shape(cdata)))}

        if block.pointData:
            for kp, pdata in block.pointData.items():
                    pinfo = {**pinfo, kp : (pdata.dtype, len(np.shape(pdata)))}

        if self.pvtktype == "pvts" :
            start = block.getIndexMin()
            end = block.getIndexMax()
        else:
            start = None
            end = None

        self.addSourceInfo(bijk=bijk,
                            bfile={bn: block.filename},
                            cinfo=cinfo,
                            pinfo=pinfo,
                            start=start,
                            end=end)


    def addSourceInfo(self, bijk, bfile, cinfo, pinfo, start=None, end=None):
        """
        Add block cell and point infos (arrays datatype), start and end indices of the blocks

        Parameters
        -----------
            bijk : array-like of int tuple
                Block number Ids
            bfile : dict
                Format : "f{blockId}": filename
                with blockId the number Id out of the total number of blocks
            cinfo : dict
                Cell array names, datatype and number of components
            pinfo : dict
                Point array names, datatype and number of components
            starts : dict, optional
                Indices min associated to the blocks
                Required for the export in PVTS format
            ends : dict, optional
                Indices max associated to the blocks
                Required for the export in PVTS format
        """
        if not self.vtkfiles:
            self.vtkfiles = bfile
        else:
            self.vtkfiles = {**self.vtkfiles, **bfile}
        
        self.addCellInfo(cinfo=cinfo)
        self.addPointInfo(pinfo=pinfo)
        self._addBlockExtents(bijk, start, end)


    def getSources(self):
        """
        Return the files of the blocks constituting the PVTK

        Returns
        --------
            dict
                Block Id and associated rootnames
        """
        return self.vtkfiles


    def getCellInfo(self):
        """
        Return the cell data array names, types and number of components
        
        Returns
        --------
            dict
                Cell data array name associated to a tuple of the type and number of components
        """
        return self.cellInfo


    def getPointInfo(self):
        """
        Return the point data array names, types and number of components
        
        Returns
        --------
            dict
                Point data array name associated to a tuple of the type and number of components
        """
        return self.pointInfo


    def getBlocksStarts(self):
        """
        Return the blocks minimal indices

        Returns
        -------
            dict
                Block ID associated to min indices tuple
        """
        if self.starts is None:
            return {}
        else:
            return self.starts


    def getBlocksEnds(self):
        """
        Return the blocks maximal indices

        Returns
        -------
            dict
                Block ID associated to max indices tuple
        """
        if self.ends is None:
            return {}
        else:
            return self.ends
    

    def getAllBlocksIndices(self):
        """
        Return the list of all block Ids (ni, nj, nk) 

        Returns
        -------
            list of tuple of int
        """
        nb1, nb2, nb3 = self.nblocks
        b1, b2, b3 = np.meshgrid(range(nb1), range(nb2), range(nb3))
        nijk = [(ni, nj, nk) for ni, nj, nk in zip(b1.reshape(-1), b2.reshape(-1), b3.reshape(-1))]
        
        return nijk


    def export(self, gids=False, writeBlocks=True):
        """
        Write the VTK files as well as the PVTK file

        Parameters
        -----------
            gids : bool, optional
                Whether or not to add the global Ids to the files
            writeBlocks : bool, optional
                If true, all blocks are exported. If False, only the PVTK header is written
        """
        if writeBlocks:
            self.__setVtkModels()
            for vtkmodel in self.vtkmodels:
                vtkmodel.export(gids=gids)

        self.writeParallelFile()


    def writeParallelFile(self):
        """Write the PVTK file of the model
        """
        nijk = self.getAllBlocksIndices()

        vtksources = [os.path.split(self.vtkfiles[bn])[-1] for bn in range(self.nbltot)]
        cdata = (self.n, np.dtype('float32')) #extent

        if self.pvtktype == "pvts":
            starts = [self.starts[(*bid, 0)] for bid in nijk]
            ends = [self.ends[(*bid, 1)] for bid in nijk]

        else:
            starts = None
            ends = None

        vtkwriter.writeParallelVTKGrid(path=self.rootname,
                                sources=vtksources,
                                coordsData=cdata,
                                starts=starts,
                                ends=ends,
                                cellData=self.cellInfo,
                                pointData=self.pointInfo)

    def getReader(self):
        """
        Return the appropriate vtk.vtkXmlReader
        """
        if self.vtktype == "pvtu":
            return vtk.vtkXMLPUnstructuredGridReader()
        elif self.vtktype == "pvts":
            return vtk.vtkXMLPStructuredGridReader()
        else:
            print("Unrecognized Parallel file format.")
            return None


    def getData(self):
        """
        Read the data model

        Returns
        -------
            vtkStructuredGrid if self.pvtktype is "pvts" or vtkUnstructuredGrid if "pvtu"
        """
        reader = self.getReader()
        reader.SetFilename(self.filename)
        reader.Update()
        
        return reader.GetOutput()
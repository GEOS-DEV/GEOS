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
import numba
from pyevtk.vtk import VtkFile, VtkUnstructuredGrid, VtkStructuredGrid, VtkPUnstructuredGrid, VtkPStructuredGrid, VtkParallelFile
from pyevtk.hl import _appendDataToFile


def _addDataToFile(vtkFile, cellData, pointData):
    """
    Modified from pyevtk.hl, Copyright 2010 - 2016 Paulo A. Herrera. All rights reserved."""
    if pointData:
        keys = list(pointData.keys())
        if None in keys:
            raise ValueError("Please check that you have correctly provided a key name for all datasets")
        # find first scalar and vector data key to set it as attribute
        gpids = next((key for key in keys if "GlobalPointIds" in key), None)
        scalars = next(
            (key for key in keys if isinstance(pointData[key] , np.ndarray) and key!=gpids), None
        )
        vectors = next((key for key in keys if isinstance(pointData[key], tuple)), None)
        
        vtkFile.openData("Point", scalars=scalars, vectors=vectors)
        if gpids:
            vtkFile.xml.addAttributes (GlobalIds=gpids)

        for key in keys:
            data = pointData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Point")

    # Cell data
    if cellData:
        keys = list(cellData.keys())
        if None in keys:
            raise ValueError("Please check that you have correctly provided a key name for all datasets")
        # find first scalar and vector data key to set it as attribute
        gcids = next((key for key in keys if ("GlobalCellIds" in key)), None)
        scalars = next(
            (key for key in keys if isinstance(cellData[key] , np.ndarray) and key!=gcids), None
        )
        vectors = next((key for key in keys if isinstance(cellData[key], tuple)), None)

        vtkFile.openData("Cell", scalars=scalars, vectors=vectors)
        if gcids:
            vtkFile.xml.addAttributes (GlobalIds=gcids)
        for key in keys:
            data = cellData[key]
            vtkFile.addData(key, data)
        vtkFile.closeData("Cell")   


def structuredToVTK(path, x, y, z, cellData = None, pointData = None, start=(0,0,0)):
    """create the .vts file
    Modified from : evtk.hl, Copyright (c) 2021 Paulo A. Herrera

    Args:
        path (str): path and name of the .vts file (without the extension)
        x,y,z (3 np.array) : arrays of point coordinates along the x, y, and z axes
        pointData (dict, optional): dictionary containing arrays with node centered data.
                Keys should be the names of the data arrays. Arrays must have same dimension
                in each direction and they should be equal to the dimensions of the cell data 
                plus one and must contain only scalar data. Default to None.
        cellData (dict, optional): dictionary containing arrays with cell centered data.
                Keys should be the names of the data arrays. Arrays must have the same dimensions 
                in all directions and must contain only scalar data. Default to None.

    Returns:
        str: Full path to saved file.
    """
    assert (x.ndim == 3 and y.ndim == 3 and z.ndim == 3), "Wrong arrays dimensions"
    
    ftype = VtkStructuredGrid
    s = x.shape
    nx, ny, nz = s[0] - 1, s[1] - 1, s[2] - 1

    # Write extent
    end = (start[0] + nx, start[1] + ny, start[2] + nz)

    w =  VtkFile(path, ftype)
    w.openGrid(start = start, end = end)
    w.openPiece(start = start, end = end)
    w.openElement("Points")
    w.addData("points", (x,y,z))
    w.closeElement("Points")

    _addDataToFile(w, pointData=pointData, cellData=cellData)
    
    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    
    _appendDataToFile(w, pointData=pointData, cellData=cellData)

    w.save()
    return w.getFileName()


def unstructuredGridToVTK(path, x, y, z, connectivity, offsets, cell_types, pointData=None, cellData=None):
    """
    Modified from pyevtk.hl, Copyright 2010 - 2016 Paulo A. Herrera. All rights reserved.

    Export unstructured grid and associated data.

    Parameters
    ----------
    path : str
        name of the file without extension where data should be saved.
    x : array-like
        x coordinates of the vertices.
    y : array-like
        y coordinates of the vertices.
    z : array-like
        z coordinates of the vertices.
    connectivity : array-like
        1D array that defines the vertices associated to each element.
        Together with offset define the connectivity or topology of the grid.
        It is assumed that vertices in an element are listed consecutively.
    offsets : array-like
        1D array with the index of the last vertex of each element
        in the connectivity array.
        It should have length nelem,
        where nelem is the number of cells or elements in the grid..
    cell_types : TYPE
        1D array with an integer that defines the cell type of
        each element in the grid.
        It should have size nelem.
        This should be assigned from evtk.vtk.VtkXXXX.tid, where XXXX represent
        the type of cell.
        Please check the VTK file format specification for allowed cell types.
    cellData : dict, optional
        dictionary with variables associated to each cell.
        Keys should be the names of the variable stored in each array.
        All arrays must have the same number of elements.
    pointData : dict, optional
        dictionary with variables associated to each vertex.
        Keys should be the names of the variable stored in each array.
        All arrays must have the same number of elements.

    Returns
    -------
    str
        Full path to saved file.
    """
    assert x.size == y.size == z.size
    x = np.array(x)
    y = np.array(y)
    z = np.array(z) 
    connectivity = np.array(connectivity)
    offsets = np.array(offsets)
    cell_types = np.array(cell_types)
    
    npoints = x.size
    ncells = cell_types.size
    assert offsets.size == ncells
    
    w = VtkFile(path, VtkUnstructuredGrid)
    w.openGrid()
    w.openPiece(ncells=ncells, npoints=npoints)
    
    w.openElement("Points")
    w.addData("points", (x, y, z))
    w.closeElement("Points")
    w.openElement("Cells")
    w.addData("connectivity", connectivity)
    w.addData("offsets", offsets)
    w.addData("types", cell_types)
    w.closeElement("Cells")

    
    _addDataToFile(w, cellData=cellData, pointData=pointData)

    w.closePiece()
    w.closeGrid()
    w.appendData( (x,y,z) )
    w.appendData(connectivity).appendData(offsets).appendData(cell_types)

    _appendDataToFile(w, pointData=pointData, cellData=cellData)

    w.save()
    return w.getFileName()


def writeParallelVTKGrid(path, sources, coordsData, starts=None, ends=None, ghostlevel=0, cellData=None, pointData=None):
    """
    Modified from pyevtk.hl, Copyright 2010 - 2016 Paulo A. Herrera. All rights reserved. 

    Writes a parallel vtk file from grid-like data:
    VTKStructuredGrid or VTKUnstructuredGrid

    Parameters
    ----------
        path : str
            name of the file without extension.
        coordsData : tuple
            2-tuple (shape, dtype) where shape is the
            shape of the coordinates of the full mesh
            and dtype is the dtype of the coordinates.
        starts : list
            list of 3-tuple representing where each source file starts
            in each dimension
        sources : list
            list of the relative paths of the source files where the actual data is found
        ghostlevel : int, optional
            Number of ghost-levels by which
            the extents in the individual source files overlap.
        pointData : dict
            dictionnary containing the information about the arrays
            containing node centered data.
            Keys shoud be the names of the arrays.
            Values are (dtype, number of components)
        cellData :
            dictionnary containing the information about the arrays
            containing cell centered data.
            Keys shoud be the names of the arrays.
            Values are (dtype, number of components)
    """
    common_ext = sources[0].split(".")[-1]
    assert all(s.split(".")[-1] == common_ext for s in sources)

    if common_ext == "vts":
        assert len(starts) == len(ends) == len(sources)
        ftype = VtkPStructuredGrid
    elif common_ext == "vtu":
        ftype = VtkPUnstructuredGrid
    else:
        raise TypeError("This function only works with VTU or VTS")


    w = VtkParallelFile(path, ftype)
    
    if common_ext == "vts":
        start = (0, 0, 0)
        (s_x, s_y, s_z), dtype = coordsData
        end = s_x - 1, s_y - 1, s_z - 1
        w.openGrid(start=start, end=end, ghostlevel=ghostlevel)
        
    elif common_ext == "vtu":
        _, dtype = coordsData
        w.openGrid(ghostlevel=ghostlevel)
    
    _addDataToParallelFile(w, cellData, pointData)

    w.openElement("PPoints")
    w.addHeader("points", dtype=dtype, ncomp=3)
    w.closeElement("PPoints")

    if common_ext == "vts":
        for start_source, end_source, source in zip(starts, ends, sources):
            w.addPiece(start_source, end_source, source)
    elif common_ext == "vtu":
        for source in sources:
            w.addPiece(source=source)

    w.closeGrid()
    w.save()
    return w.getFileName()

def _addDataToParallelFile(vtkParallelFile, cellData, pointData):
    """
    Modified from pyevtk.hl, Copyright 2010 - 2016 Paulo A. Herrera. All rights reserved.
    """
    assert isinstance(vtkParallelFile, VtkParallelFile)
    # Point data
    if pointData:
        keys = list(pointData.keys())
        # find first scalar and vector data key to set it as attribute
        gpids = next((key for key in keys if "GlobalPointIds" in key), None)
        scalars = next((key for key in keys if pointData[key][1] == 1), None)
        vectors = next((key for key in keys if pointData[key][1] == 3), None)
    
        vtkParallelFile.openData("PPoint", scalars=scalars, vectors=vectors)
        if gpids:
            vtkParallelFile.xml.addAttributes (GlobalIds=gpids)
        for key in keys:
            dtype, ncomp = pointData[key]
            vtkParallelFile.addHeader(key, dtype=dtype, ncomp=ncomp)
        vtkParallelFile.closeData("PPoint")

    # Cell data
    if cellData:
        keys = list(cellData.keys())
        # find first scalar and vector data key to set it as attribute
        gcids = next((key for key in keys if "GlobalCellIds" in key), None)
        scalars = next((key for key in keys if cellData[key][1] == 1), None)
        vectors = next((key for key in keys if cellData[key][1] == 3), None)
        
        vtkParallelFile.openData("PCell", scalars=scalars, vectors=vectors)
        if gcids:
            vtkParallelFile.xml.addAttributes (GlobalIds=gcids)
        for key in keys:
            dtype, ncomp = cellData[key]
            vtkParallelFile.addHeader(key, dtype=dtype, ncomp=ncomp)
        vtkParallelFile.closeData("PCell")




@numba.jit(nopython=True)
def xyz(Nx,Ny,Nz,dx,dy,dz,xmin,ymin,zmin):
    """ Compute points coordinates (vertices) for .vtu file

    Args:
        Nx (int): number of points on the x axis
        Ny (int): number of points on the y axis
        Nz (int): number of points on the z axis
        dx (float): x-axis spacing
        dy (float): y-axis spacing
        dz (float): z-axis spacing
        xmin (float): minimum x value
        ymin (float): minimum y value
        zmin (float): minimum z value

    Returns:
        x,y,z (3 np.array) : arrays of point coordinates along the x, y, and z axes
    """
    x=np.zeros(Nz*Nx*Ny, dtype='float32')
    y=np.zeros(Nz*Nx*Ny, dtype='float32')
    z=np.zeros(Nz*Nx*Ny, dtype='float32')
    for k in range(0,Nz):
        k_index=k*dz+zmin
        for j in range(0,Ny):
            j_index=j*dy+ymin
            for i in range(0,Nx):
                index=i+j*Nx+k*Ny*Nx
                x[index]=i*dx+xmin
                y[index]=j_index
                z[index]=k_index
    return x,y,z


@numba.jit(nopython=True)
def x_y_z(nx,ny,nz,dx,dy,dz,xmin,ymin,zmin):
    """ Compute points coordinates (vertices) for .vts file

    Args:
        nx (int): number of points on the x axis
        ny (int): number of points on the y axis
        nz (int): number of points on the z axis
        dx (float): x-axis spacing
        dy (float): y-axis spacing
        dz (float): z-axis spacing
        xmin (float): minimum x value
        ymin (float): minimum y value
        zmin (float): minimum z value

    Returns:
        x,y,z (3 np.array) : arrays of point coordinates along the x, y, and z axes
    """
    x = np.zeros((nx, ny, nz), dtype='float32')
    y = np.zeros((nx, ny, nz), dtype='float32')
    z = np.zeros((nx, ny, nz), dtype='float32')
    for k in range(0,nz):
        k_index=k*dz+zmin
        for j in range(0,ny):
            j_index=j*dy+ymin
            for i in range(0,nx):
                x[i,j,k]=i*dx+xmin
                y[i,j,k]=j_index
                z[i,j,k]=k_index
    return x,y,z


@numba.jit(nopython=True)
def connectivity(Nx,Ny,Nz):
    """Compute connectivity and offsets

    Args:
        Nx (int): number of points on the x axis
        Ny (int): number of points on the y axis
        Nz (int): number of points on the z axis

    Returns:
        conn, offset (2 np.arrays):  2 arrays of connectivity and offsets
        numCells (int): number of cells
    """
    numCells=(Nx-1)*(Ny-1)*(Nz-1)
    conn = np.zeros(numCells*8, dtype='int64')
    offset=np.zeros(numCells, dtype='int64')
    numPoint=np.zeros(8, dtype='uint32')
    cellIndex=0
    for k in range(0,Nz-1):
        for j in range(0,Ny-1):
            for i in range(0,Nx-1):
                offsetIndex=8*(i+j*(Nx-1)+k*(Ny-1)*(Nx-1))
                numPoint[0]=i+j*Nx+k*Ny*Nx
                numPoint[1]=i+1+j*Nx+k*Ny*Nx
                numPoint[2]=i+1+(j+1)*Nx+k*Ny*Nx
                numPoint[3]=i+(j+1)*Nx+k*Ny*Nx
                numPoint[4]=i+j*Nx+(k+1)*Ny*Nx
                numPoint[5]=i+1+j*Nx+(k+1)*Ny*Nx
                numPoint[6]=i+1+(j+1)*Nx+(k+1)*Ny*Nx
                numPoint[7]=i+(j+1)*Nx+(k+1)*Ny*Nx
                conn[offsetIndex]=numPoint[0]
                conn[offsetIndex+1]=numPoint[1]
                conn[offsetIndex+2]=numPoint[2]
                conn[offsetIndex+3]=numPoint[3]
                conn[offsetIndex+4]=numPoint[4]
                conn[offsetIndex+5]=numPoint[5]
                conn[offsetIndex+6]=numPoint[6]
                conn[offsetIndex+7]=numPoint[7]
                offset[cellIndex]=int(offsetIndex+8)
                cellIndex+=1
    return conn,offset,numCells

@numba.jit(nopython=True)
def pGlobalIds(Nx,Ny,Nz,dx,dy,dz,xmin,ymin,zmin,glNx,glNy,glNz,xmin0,ymin0,zmin0):
    """Compute global ids for points

    Args:
        Nx (int): number of points on the x axis for the block
        Ny (int): number of points on the y axis for the block
        Nz (int): number of points on the z axis for the block
        dx (float): x-axis spacing
        dy (float): y-axis spacing
        dz (float): z-axis spacing
        xmin (float): minimum x value for the block
        ymin (float): minimum y value for the block
        zmin (float): minimum z value for the block
        glNx (int): global number of points on the x axis
        glNy (int): global number of points on the y axis
        glNz (int): global number of points on the z axis
        xmin0 (float): minimum x value for the _0 file
        ymin0 (float): minimum y value for the _0 file
        zmin0 (float): minimum z value for the _0 file

    Returns:
        pIds (numpy.array): point global ids for the block
    """
    offx=int(round((xmin-xmin0)/dx))
    offy=int(round((ymin-ymin0)/dy))
    offz=int(round((zmin-zmin0)/dz))
    pIds=np.zeros(Nz*Nx*Ny, dtype='int64')
    for k in range(0,Nz):
        for j in range(0,Ny):
            for i in range(0,Nx):
                index=i+j*Nx+k*Ny*Nx
                glindex=(i+offx)+(j+offy)*glNx+(k+offz)*glNy*glNx
                pIds[index]=glindex
    return pIds

@numba.jit(nopython=True)
def cGlobalIds(Nx,Ny,Nz,dx,dy,dz,xmin,ymin,zmin,glNx,glNy,glNz,xmin0,ymin0,zmin0):
    """Compute global ids for cells

    Args:
        Nx (int): number of points on the x axis for the block
        Ny (int): number of points on the y axis for the block
        Nz (int): number of points on the z axis for the block
        dx (float): x-axis spacing
        dy (float): y-axis spacing
        dz (float): z-axis spacing
        xmin (float): minimum x value for the block
        ymin (float): minimum y value for the block
        zmin (float): minimum z value for the block
        glNx (int): global number of points on the x axis
        glNy (int): global number of points on the y axis
        glNz (int): global number of points on the z axis
        xmin0 (float): minimum x value for the _0 file
        ymin0 (float): minimum y value for the _0 file
        zmin0 (float): minimum z value for the _0 file

    Returns:
        cIds (numpy.array): cell global ids for the block
    """
    offx=int(round((xmin-xmin0)/dx))
    offy=int(round((ymin-ymin0)/dy))
    offz=int(round((zmin-zmin0)/dz))
    cIds=np.zeros((Nz-1)*(Nx-1)*(Ny-1), dtype='int64')
    for k in range(0,Nz-1):
        for j in range(0,Ny-1):
            for i in range(0,Nx-1):
                cellIndex=(i+j*(Nx-1)+k*(Ny-1)*(Nx-1))
                cellGlIndex=((i+offx)+(j+offy)*(glNx-1)+(k+offz)*(glNy-1)*(glNx-1))
                cIds[cellIndex]=cellGlIndex
    return cIds
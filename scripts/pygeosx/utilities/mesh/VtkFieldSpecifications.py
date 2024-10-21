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

import copy
import numpy as np
from vtk.util import numpy_support as VN


class VTKFieldSpecifications:
    """
    Object containing field dataset of a VTK format mesh

    Attributes
    ----------
        arrays : dict
            Dict containing the {array names : array} of the given dataset 
    """
    def __init__(self, fieldData):
        """
        Parameters
        ----------
            fieldData : vtk.vtkFieldData
               Data contained in the VTK mesh
        """
        self.arrays = {}
        assert(fieldData.IsA("vtkFieldData"))
        for i in range(fieldData.GetNumberOfArrays()):
            value = fieldData.GetArray(i)
            key = value.GetName()
            self.arrays.update({key: value})
        self.fieldtype = None


    def hasArray(self, name):
        """
        Check if the object contains an array with the requested name
        
        Parameters
        ----------
            name : str
                Name of the cell/point data array
        
        Returns
        -------
            bool : True if the array exists. False otherwise
        """
        if name in self.arrays.keys():
            return True
        else:
            return False


    def getArray(self, name, sorted=False):
        """ 
        Return the requested array from the cell data
        
        Parameters
        ----------
            name : str
                Name of the cell/point data array
            sorted : bool
                Sort with global ids \
                Require the presence of a global ids array
        
        Returns
        -------
            numpy array : requested array
        """
        array = None

        if self.hasArray(name):
            array = VN.vtk_to_numpy(self.arrays[name])
             
            if sorted:
                ftype = self.fieldtype.split("Data")[0]
                if self.hasArray(f"Global{ftype}Ids"):
                    gids = self.getCopyArray(f"Global{ftype}Ids")
                    array = array[np.argsort(gids)]

        return array


    def getCopyArray(self, name, **kwargs):
        """
        Return a copy of the requested array from the cell data
        
        Parameters
        ----------
            name : str
                Name of the cell/point data array
        
        Returns
        -------
            numpy array : copy of the requested array
        """

        array = self.getArray(name, **kwargs)

        if array is not None:
            array = copy.deepcopy(array)

        return array


    def getVtkArray(self, name):
        """
        Return the vtkDataArray requested
        
        Parameters
        ----------
            name : str
                Name of the cell/point data array
        
        Returns
        -------
            numpy array : copy of the requested array
        """
        if self.hasArray(name):
            return self.arrays[name]
        else:
            return None


    def setArray(self, name, value, overwrite=False):
        """
        Return a copy of the requested array from the cell data
        
        Parameters
        ----------
            name : str
                Name of the cell data array
        
        Returns
        -------
            numpy array : copy of the requested array
        """
        if self.hasArray(name) and overwrite == False:
            print(f"Warning! \n This dataset already contains a cell data array named {name}. Set the 'overwrite' parameter to True to bypass this warning")
        else:
            array = VN.vtk_to_numpy(self.arrays[name])
            array[:] = value[:]



class VTKCellSpecifications(VTKFieldSpecifications):
    """
    Contains the cell data information from a VTK Mesh
    Inherits from VTKFieldSpecifications
    """
    def __init__(self, celldata):
        """
        Parameters
        ----------
            celldata : vtk.vtkCellData
                Cell data of the mesh
        """
        assert(celldata.IsA("vtkCellData"))
        super().__init__(fieldData=celldata)
        self.fieldtype = "CellData"


class VTKPointSpecifications(VTKFieldSpecifications):
    """ 
    Contains the point data information from a VTK Mesh
    Inherits from VTKFieldSpecifications

    Parameters
    ----------
        pointdata : vtk.vtkPointData
            Point data of the mesh

    Attributes
    ---------
        arrays : dict
            Dict containing the {name, vtkDataArray} of each point data array    
    """
    def __init__(self, pointdata):
        """
        Parameters
        ----------
            pointdata : vtk.vtkPointData
                Point data of the mesh
        """
        assert(pointdata.IsA("vtkPointData"))
        super().__init__(fieldData=pointdata)
        self.fieldtype = "PointData"

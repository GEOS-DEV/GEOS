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
import ctypes as ct

import pygeosx

from .WaveSolver import WaveSolver


class ElasticSolver(WaveSolver):
    """
    ElasticSolver Object containing all methods to run ElasticSEM simulation with GEOSX

    Attributes
    -----------
        dt : float
            Time step for simulation
        minTime : float
            Min time to consider
        maxTime : float
            End time to consider
        dtSeismo : float
            Time step to save pressure for seismic trace
        minTimeSim : float
            Starting time of simulation
        maxTimeSim : float
            End Time of simulation
        dtWaveField : float
            Time step to save fields
        sourceType : str
            Type of source
        sourceFreq : float
            Frequency of the source
        name : str
            Solver name
        type : str
            Type of solver
            Default is None
        geosxArgs : GeosxArgs
            Object containing GEOSX launching options
        geosx : pygeosx instance
            
        alreadyInitialized : bool
            Flag indicating if the initial conditions have been applied yet
        firstGeosxInit : bool
            Flag for initialize or reinitialize 
        collectionTargets : list
            Output targets for geosx
        hdf5Targets : list
            HDF5 output targets for geosx
        vtkTargets : list
            VTK output targets for geosx
    """
    def __init__(self,
                 dt=None,
                 minTime=0,
                 maxTime=None,
                 dtSeismo=None,
                 dtWaveField=None,
                 sourceType=None,
                 sourceFreq=None,
                 **kwargs):
        """
        Parameters
        ----------
            dt : float
                Time step for simulation
            minTime : float
                Starting time of simulation
                Default is 0
            maxTime : float
                End Time of simulation
            dtSeismo : float
                Time step to save pressure for seismic trace
            dtWaveField : float
                Time step to save fields
            sourceType : str
                Type of source
                Default is None
            sourceFreq : float
                Frequency of the source
                Default is None
            kwargs : keyword args
                geosx_argv : list
                    GEOSX arguments or command line as a splitted line
        """
        super().__init__(dt=dt,
                         minTime=minTime,
                         maxTime=maxTime,
                         dtSeismo=dtSeismo,
                         dtWaveField=dtWaveField,
                         sourceType=sourceType,
                         sourceFreq=sourceFreq,
                         **kwargs)
        
    def initialize(self, rank=0, xml=None):
        super().initialize(rank, xml)
        try:
            useDAS = self.xml.getAttribute(parentElement=self._getType(), attributeTag="useDAS")

        except AttributeError:
            useDAS = None
        
        if useDAS == "none":
            try:
                linearGEO = bool(self.xml.getAttribute(self._getType(), "linearDASGeometry"))
            except AttributeError:
                linearGEO = False
        
            if linearGEO is True:
                self.useDAS = True


    def __repr__(self):
        string_list = []
        string_list.append("Solver type : " + type(self).__name__ + "\n")
        string_list.append("dt : " + str(self.dt) + "\n")
        string_list.append("maxTime : " + str(self.maxTime) + "\n")
        string_list.append("dtSeismo : " + str(self.dtSeismo) + "\n")
        string_list.append("Outputs : " + str(self.hdf5Targets) + "\n" + str(self.vtkTargets) + "\n")
        rep = ""
        for string in string_list:
            rep += string

        return rep


    def updateVelocityModel(self, vel, component, **kwargs):
        """
        Update velocity value in GEOS

        Parameters
        ----------
            vel : float/array
                Value(s) for velocity field
            component : str
                Vs or Vp
        """
        assert component.lower() in ("vs", "vp"), "Only Vs or Vp component accepted"
        super().updateVelocityModel(vel, velocityName="elasticVelocity"+component.title(), **kwargs)


    def getVelocityModel(self, component, filterGhost=False, **kwargs):
        """
        Get the velocity values

        Parameters
        -----------
            component : str
                Vs or Vp
            filterGhost : bool
                Filter the ghost ranks

        Returns
        -------
            velocity : numpy array
                Array containing the velocity values
        """
        assert component.lower() in ("vs", "vp"), "Only Vs or Vp component accepted"
        
        velocity = super().getVelocityModel(velocityName="elasticVelocity"+component.title(), filterGhost=filterGhost, **kwargs)

        return velocity


    def getDensityModel(self, filterGhost=False, **kwargs):
        """
        Get the density values

        Parameters
        -----------
            filterGhost : bool
                Filter the ghost ranks

        Returns
        --------
            density : numpy array
                Array containing the density values
        """
        density = self.getField("elasticDensity", filterGhost=filterGhost, **kwargs)

        return density


    def updateDensityModel(self, density, **kwargs):
        """
        Update density values in GEOS
        
        Parameters
        -----------
            density : array
                New values for the density
        """
        super().updateDensityModel( density=density, densityName="elasticDensity", **kwargs )


    def getDASSignalAtReceivers(self):
        """
        Get the DAS signal values at receivers coordinates

        Returns
        --------
            dassignal : numpy array
                Array containing the DAS signal values at all time step at all receivers coordinates
        """
        if self.type != "ElasticSEM":
          raise TypeError(f"DAS signal not implemented for solver of type {self.type}.")
        else:
          dassignal = self.solver.get_wrapper(f"dasSignalNp1AtReceivers").value().to_numpy()

        return dassignal


    def getDisplacementAtReceivers(self, component="X"):
        """
        Get the displacement values at receivers coordinates for a given direction

        Returns
        --------
            displacement : numpy array
                Array containing the displacements values at all time step at all receivers coordinates
        """
        assert component.upper() in ("X", "Y", "Z")
        if self.type == "ElasticFirstOrderSEM":
            displacement = self.solver.get_wrapper(f"displacement{component.lower()}Np1AtReceivers").value().to_numpy()
        elif self.type == "ElasticSEM":
            displacement = self.solver.get_wrapper(f"displacement{component.upper()}Np1AtReceivers").value().to_numpy()

        return displacement
    

    def getAllDisplacementAtReceivers(self):
        """
        Get the displacement for the x, y and z directions at all time step and all receivers coordinates

        Returns
        --------
            displacementX : numpy array
                Component X of the displacement
            displacementY : numpy array
                Component Y of the displacement
            displacementZ : numpy array
                Component Z of the displacement
        """
        displacementX = self.getDisplacementAtReceivers("X")
        displacementY = self.getDisplacementAtReceivers("Y")
        displacementZ = self.getDisplacementAtReceivers("Z")

        return displacementX, displacementY, displacementZ


    def resetWaveField(self, **kwargs):
        """Reinitialize all displacement values on the Wavefield to zero in GEOSX"""

        self.geosx.get_wrapper("Solvers/"+self.name+"/indexSeismoTrace").value()[0] = 0
        meshName = self._getMeshName()
        discretization = self._getDiscretization()

        nodeManagerPath = f"domain/MeshBodies/{meshName}/meshLevels/{discretization}/nodeManager/"

        if self.type == "ElasticSEM":
            for component in ("x", "y", "z"):
                for ts in ("nm1", "n", "np1"):
                    displacement = self.geosx.get_wrapper(nodeManagerPath+f"displacement{component}_{ts}").value()
                    displacement.set_access_level(pygeosx.pylvarray.MODIFIABLE)

                    displacement.to_numpy()[:] = 0.0

        elif self.type == "ElasticFirstOrderSEM":
            component = ("x", "y", "z")
            for c in component:
                displacement_np1 = self.geosx.get_wrapper(nodeManagerPath+f"displacement{c}_np1").value()
                displacement_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

                displacement_np1.to_numpy()[:] = 0.0
            
            prefix = self._getPrefixPath(**kwargs)
            for i, c in enumerate(component):
                for j in range(i, len(component)):
                    cc = c + component[j]

                sigma = self.solver.get_wrapper(prefix+f"stresstensor{cc}").value()
                sigma.set_access_level(pygeosx.pylvarray.MODIFIABLE)

                sigma.to_numpy()[:] = 0.0


    def resetDisplacementAtReceivers(self):
        """Reinitialize displacement values at receivers to 0
        """
        for component in ("X", "Y", "Z"):
            displacement = self.solver.get_wrapper(f"displacement{component}Np1AtReceivers").value()
            displacement.set_access_level(pygeosx.pylvarray.MODIFIABLE)
 
            displacement.to_numpy()[:] = 0.0


    def getWaveField( self ):
        if self.useDAS:
            return self.getDASSignalAtReceivers()
        else:
            return self.getAllDisplacementAtReceivers()


    def getFullWaveFieldAtReceivers(self, comm):
        print( "This method is not implemented yet" )

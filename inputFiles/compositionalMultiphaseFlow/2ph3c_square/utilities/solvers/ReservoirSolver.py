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

from .Solver import Solver


class ReservoirSolver(Solver):
    """
    Reservoir solver object containing methods to run reservoir simulations in GEOS

    Attributes
    -----------
        xml : XML
            XML object containing parameters for GEOSX initialization
        initialDt : float
            Time step for the simulation
        maxTime : float
            End time of simulation
        name : str
            Solver name
        type : str
            Type of solver
            Default is None
        geosxArgs : GeosxArgs
            Object containing GEOSX launching options
        geosx : pygeosx.Group
            pygeosx initialize instance
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
                initialDt=None,
                maxTime=None,
                **kwargs):
        """
        Parameters
        -----------
            initialDt : float
                Initial time step \
                Default is None
            maxTime : float
                End time of simulation \
                Default is None
        """
        super().__init__(**kwargs)

        self.initialDt = initialDt
        self.maxTime = maxTime

        self.isCoupled = kwargs.get("coupled", False)
        self.dtSeismic4D = None


    def _setTimeVariables(self):
        """Set the time variables attributes"""
        # Set up variables from xml
        events = self.xml.events
        outputs = self.xml.outputs

        if self.isCoupled:
            if self.dtSeismic4D is None:
                for event in events["PeriodicEvent"]:
                    if event["name"] == "seismic4D":
                        self.dtSeismic4D = float(event["timeFrequency"])
            else:
                self.updateTimeVariable("dtSeismic4D")
    
        if self.maxTime is None:
            self.maxTime = float(events["maxTime"])
        else:
            self.updateTimeVariable("maxTime")

        fvm = self.xml.solvers[self.type] 
        if self.initialDt is None:
            for k, v in fvm.items():
                if k == "initialDt":
                    self.initialDt = np.array(v, dtype=float)
        else:
            self.updateTimeVariable("initialDt")

  
    def updateTimeVariable(self, variable):
        """Overwrite a time variable in GEOS"""
        if variable == "maxTime":
            assert hasattr(self, "maxTime")
            self.geosx.get_wrapper("Events/maxTime").value()[0] = self.maxTime
        elif variable == "initialDt":
            assert hasattr(self, "initialDt")
            self.geosx.get_wrapper(f"/Solvers/{self.name}/initialDt").value()[0] = self.initialDt
        elif variable == "dtSeismic4D":
            assert hasattr(self, "dtSeismic4D")
            self.geosx.get_wrapper("Events/seismic4D/timeFrequency").value()[0] = self.dtSeismic4D
	

    def execute(self, time):
        """
        Do one solver iteration

        Parameters
        ----------
            time : float
                Current time of simulation
            dt : float
                Timestep
        """
        self.solver.execute(time, self.initialDt)


    def getTimeStep(self):
        """
        Get the value of `initialDt` variable from GEOS

        Returns
        --------
            float
                initialDt value
        """
        return self.solver.get_wrapper("initialDt").value()[0]
    

    def updateTimeStep(self):
        """Update the attribute value with the one from GEOS"""
        self.initialDt = self.getTimeStep()


    def getPressure(self, **kwargs):
        """
        Get the local pressure
        
        Returns
        --------
            pressure : numpy array
                Local pressure
        """
        pressure = self.getField("pressure", **kwargs)
    
        return pressure


    def getDeltaPressure(self, **kwargs):
        """
        Get the local delta pressure
        
        Returns
        --------
            deltaPressure : numpy array
                Local delta pressure
        """
        deltaPressure = self.getField("deltaPressure", **kwargs)
    
        return deltaPressure


    def getPhaseVolumeFractionGas(self, **kwargs):
        """
        Get the local gas phase volume fraction
        
        Returns
        --------
            phaseVolumeFractionGas : numpy array
                Local gas phase volume fraction
        """
        phaseVolumeFraction = self.getField("phaseVolumeFraction", **kwargs)
        phaseVolumeFractionGas = np.ascontiguousarray(phaseVolumeFraction[:,0])

        return phaseVolumeFractionGas


    def getPhaseVolumeFractionWater(self, **kwargs):
        """
        Get the local water phase volume fraction
        
        Returns
        --------
            phaseVolumeFractionWater : numpy array
                Local water phase volume fraction
        """
        phaseVolumeFraction = self.getField("phaseVolumeFraction", **kwargs)
        phaseVolumeFractionWater = np.ascontiguousarray(phaseVolumeFraction[:,1])

        return phaseVolumeFractionWater


    def getRockPorosity(self, **kwargs):
        """
        Get the local rock porosity
        
        Returns
        --------
            rockPorosity : numpy array
                Local rock porosity
        """
        rockPorosity = self.getField("rockPorosity_referencePorosity", **kwargs)

        return rockPorosity

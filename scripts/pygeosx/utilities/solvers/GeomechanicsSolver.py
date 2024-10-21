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

from .Solver import Solver


class GeomechanicsSolver(Solver):
    """
    Geomechanics solver object containing methods to run geomechanics simulations in GEOS

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
                dt=None,
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

        self.dt = dt
        self.maxTime = maxTime


    def _setTimeVariables(self):
        """Set the time variables attributes"""
        # Set up variables from xml
        events = self.xml.events

        if self.maxTime is None:
            self.maxTime = float(events["maxTime"])
        else:
            self.updateTimeVariable("maxTime")

        if self.dt is None:
            for event in events["PeriodicEvent"]:
                if isinstance(event, dict):
                    poroName = "poromechanicsSolver"
                    if event["target"] == "/Solvers/" + poroName:
                        self.dt = float(event['forceDt'])
                else:
                    if event == "target" and events["PeriodicEvent"]["target"] == "/Solvers/" + self.name:
                        self.dt = float(events["PeriodicEvent"]["forceDt"])
        else:
                self.updateTimeVariable("dt")

  
    def updateTimeVariable(self, variable):
        """Overwrite a time variable in GEOS"""
        if variable == "maxTime":
            assert hasattr(self, "maxTime")
            self.geosx.get_wrapper("Events/maxTime").value()[0] = self.maxTime
        elif variable == "dt":
            assert hasattr(self, "dt")
            self.geosx.get_wrapper("Events/solverApplications/forceDt").value()[0]   


    def initialize( self, rank=0, xml=None ):
        super().initialize( rank, xml, stype="SinglePhasePoromechanics" )

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

import pygeosx
from ..input.Xml import XML
from ..input.GeosxArgs import GeosxArgs

from mpi4py import MPI

class Solver:
    """
    Solver class containing the main methods of a GEOS solver
    """

    def __init__(self, **kwargs):
        self.alreadyInitialized = False
        self.type = None

        argv = kwargs.get("geosx_argv", sys.argv)
        self.geosxArgs = GeosxArgs(argv)
        self.xml = None
   
 
    def initialize(self, rank=0, xml=None, **kwargs):
        """Initialization or reinitialization of GEOSX

        Parameters
        ----------
            rank : int
                Process rank
            xml : XML
                XML object containing parameters for GEOSX initialization.
                Only required if not set in the __init__ OR if different from it
        """
        if xml:
            self.updateXml(xml)
        else:
            if self.xml is None:
                try:
                    self.xml = XML(self.geosxArgs.options["xml"])
                except:
                    raise ValueError("You need to provide a xml input file")


        if not self.alreadyInitialized:
            # if GEOS is UNINITIALIZED
            if self.getGEOSState() == 0:
                self.geosx = pygeosx.initialize(rank, self.geosxArgs.getCommandLine())
                self.alreadyInitialized = True

            # elif GEOS is INITIALIZED OR READY TO RUN
            elif self.getGEOSState() in (1, 2):
                self.geosx = pygeosx.reinit(self.geosxArgs.getCommandLine())
                self.alreadyInitialized = True

            # else if COMPLETED state
            else:
                print(f"The current state of GEOS does not allow for initialization\nCurrent state: {self.getGEOSState()}")

            self.name = self._getName()
            stype = kwargs.get( "stype", None )
            self._setSolverGroup( stype )

            # Check all time variables are set/ set them from xml
            self._setTimeVariables()

            # Set the outputs collections/targets
            self.__setOutputs()


    def _setSolverGroup( self, stype=None ):
        if stype is None:
            if self.name is not None:
                if isinstance( self.name, str ):
                    self.solver = self.geosx.get_group( "/Solvers/" + self.name )
        else:
            name = self._getName( stype=stype )
            self.solver = self.geosx.get_group( "/Solvers/" + name )


    def getGEOSState(self):
        """
        Return the current GEOS state

        Returns
        ---------
            int
                GEOS state
                0 : UNINITIALIZED
                1 : INITIALIZED
                2 : READY TO RUN
                3 : COMPLETED
        """
        return pygeosx.getState()


    def _setTimeVariables(self):
        """Initialize the time variables. Specific to each solver"""
        pass

            
    def __setOutputs(self):
        if hasattr(self.xml, "outputs"):
            outputs = self.xml.outputs
            self.__setOutputTargets(outputs)

            self.collections = []
            self.hdf5Outputs = []
            self.vtkOutputs = []

            if not outputs == None:
                # Set the collections
                for target in self.collectionTargets:
                    self.collections.append(self.geosx.get_group(target))

                for target in self.hdf5Targets:
                    self.hdf5Outputs.append(self.geosx.get_group(target))

                for target in self.vtkTargets:
                    self.vtkOutputs.append(self.geosx.get_group(target))


    def __setOutputTargets(self, outputs):
        """
        Set the list of outputs targets

        Parameters
        -----------
            outputs : dict
                Dictionnary extracted from the xml input file
        """
        self.collectionTargets = []
        self.hdf5Targets = []
        self.vtkTargets = []

        if outputs:
            if isinstance(list(outputs.values())[0], list):
                if "TimeHistory" in outputs.keys():
                    for hdf5 in outputs["TimeHistory"]:
                        self.collectionTargets.append(hdf5['sources'].strip("{} "))
                        self.hdf5Targets.append("Outputs/"+hdf5['name'])

                if "VTK" in outputs.keys():
                    for vtk in outputs["VTK"]:
                        self.vtkTargets.append("Outputs/"+vtk['name'])

            else:
                if "TimeHistory" in list(outputs.keys()):
                    hdf5 = outputs["TimeHistory"]
                    self.collectionTargets.append(hdf5['sources'].strip("{} "))
                    self.hdf5Targets.append("Outputs/"+hdf5['name'])

                if "VTK" in list(outputs.keys()):
                    vtk = outputs["VTK"]
                    self.vtkTargets.append("Outputs/"+vtk['name'])


    def _getType(self):
        """
        Set the type of solver given in the xml 'Solvers' block
        
        Raises
        -------
            ValueError : if no solver is provided in the XML
        """
        if self.xml is not None:
            typesOfSolvers = self.xml.getSolverType()

            if len( typesOfSolvers ) == 1:
                self.type = typesOfSolvers[0]
            elif len( typesOfSolvers ) > 1:
                self.type = typesOfSolvers
            else:
                raise ValueError("You must provide a Solver in the XML input file")


    def _getName(self, stype=None):
        """
        Get the solver 'name' attribute from the xml
        
        Returns
        -------
            str or list of str
                Solver name from the xml \
                List of solvers name if several solvers in xml
        """
        if self.xml is not None:
            if stype is None:
                if self.type is None:
                    self._getType()
                stype = self.type
    
            # Check only one type of solver 
            if isinstance(stype, str):
                return self.xml.solvers[stype]["name"]

            elif isinstance( stype, list ):
                return [ self.xml.solvers[solvertype]["name"] for solvertype in stype ]


    def _getMeshName(self):
        """
        Get the mesh 'name' attribute from the xml
        
        Returns
        -------
            str
                Mesh name from the xml
        """
        if self.xml is not None:
            meshes = [m for m in self.xml.mesh]
            if len(meshes) <= 1:
                return self.xml.mesh[meshes[0]]["name"]      


    def _getDiscretization(self):
        """Get discretization from the XML

        Returns
        --------
            discretization : str
                Name of the discretization method
        """
        if self.xml is not None:
            if self.type is None:
                self._getType()

            if isinstance(self.type, str):
                return self.xml.solvers[self.type]["discretization"]
            elif isinstance(self.type, list):
                return [self.xml.solvers[solverType]["discretization"] for solverType in self.type]
                

    def _getTargetRegion(self):
        """
        Get the target region name from the xml
        
        Returns
        -------
            str
                target region from the xml
        """
        if self.xml is not None:
            if self.type is None:
                self._getType()
            
            if isinstance(self.type, str):
                targetRegionRaw = self.xml.solvers[self.type]["targetRegions"]
                targetRegion = targetRegionRaw.strip("{ }").split(",")
            
                if len(targetRegion) <= 1:
                    return targetRegion[0]


    def _getCellBlock(self):
        """
        Get the cell blocks names from the xml
        
        Returns
        -------
            str
                cell blocks from the xml
        """
        if self.xml is not None:
            cellElementRegion = self.xml.elementRegions["CellElementRegion"][0]
            cellBlocks = cellElementRegion["cellBlocks"].strip("{ }").split(",")

            if len(cellBlocks) <= 1:
                return cellBlocks[0] 
 
        
    def reinitSolver(self):
        """Reinitialize Solver"""
        self.solver.reinit()


    def applyInitialConditions(self):
        """Apply the initial conditions after GEOS (re)initialization"""
        if self.getGEOSState() == 1:
            pygeosx.apply_initial_conditions()



    def finalize(self):
        """Terminate GEOSX"""
        pygeosx._finalize()

         
    def updateXml(self, xml):
        """
        Update XML

        Parameters
        -----------
            xml : XML
                XML object corresponding to GEOSX input
        """
        self.xml = xml

        if self.geosxArgs.updateArg("xml", xml.filename):
            self.alreadyInitialized = False


    def updateHdf5OutputsName(self, directory, filenames, reinit=False):
        """
        Overwrite GEOSX hdf5 Outputs paths that have been read in the XML.

        Parameters
        ----------
            list_of_output : list of str
                List of requested output paths
            reinit : bool
                Perform reinitialization or not. Must be set to True if called after applyInitialConditions()
        """

        if not len(self.hdf5Outputs):
            raise ValueError("No HDF5 Outputs specified in XML.")
        else:
            for i in range(len(filenames)):
                os.makedirs(directory, exist_ok=True)

                self.hdf5Outputs[i].setOutputName(os.path.join(directory, filenames[i]))
                if reinit:
                    self.hdf5Outputs[i].reinit()


    def updateVtkOutputsName(self, directory):
        """
        Overwrite GEOSX vtk Outputs paths that have been read in the XML.

        Parameters
        ----------
            list_of_output : list of str
                List of vtk output paths
            reinit : bool
                Perform reinitialization or not. Must be set to True if called after applyInitialConditions()
        """
        if not len(self.vtkOutputs):
            pass
        else:
            self.vtkOutputs[0].setOutputDir(directory)


    def execute(self, time):
        """
        Do one solver iteration

        Parameters
        ----------
            time : float
                Current time of simulation
        """

        self.solver.execute(time, self.dt)


    def cleanup(self, time):
        """
        Finalize simulation. Also triggers write of leftover seismogram data

        Parameters
        ----------
        time : float
            Current time of simulation
        """
        self.solver.cleanup(time)


    def outputVtk(self, time):
        """
        Trigger the VTK output
 
        Parameters
        ----------
            time : float
                Current time of simulation
        """
        for vtkOutput in self.vtkOutputs:
            vtkOutput.output(time, self.dt)


    def _getPrefixPath(self, targetRegion=None, meshName=None, cellBlock=None):
        """
        Return the prefix path to get wrappers or fields in GEOS

        Parameters
        -----------
            targetRegion : str, optional
                Name of the target Region \
                Default value is taken from the xml
            meshName : str, optional
                Name of the mesh \
                Default value is taken from the xml
            cellBlock : str, optional
                Name of the cell blocks \
                Default value is taken from the xml

        Returns
        -------
            prefix : str
                Prefix path

        Raises
        -------
            AssertionError : if the variables 'targetRegion', 'meshName' \
                or `cellBlock` have multiple or no values
        """
        if targetRegion is None:
            targetRegion = self._getTargetRegion()
        if meshName is None:
            meshName = self._getMeshName()
        if cellBlock is None:
            cellBlock = self._getCellBlock()

        discretization=self._getDiscretization()
        if discretization is None:
            discretization="Level0"
        assert None not in (targetRegion, meshName, cellBlock, discretization), "No values or multiple values found for `targetRegion`, `meshName` and `cellBlock` arguments"

        prefix = os.path.join("/domain/MeshBodies", meshName, "meshLevels", discretization, "ElementRegions/elementRegionsGroup", targetRegion, "elementSubRegions", cellBlock, "")
        return prefix
    

    def getField(self, fieldName, **kwargs):
        """
        Get the requested field as numpy array

        Parameters
        -----------
            fieldName : str
                Name of the field in GEOSX

        Returns
        -------
            field : np.array
                Field requested
        """
        prefix = self._getPrefixPath(**kwargs)
        field = self.solver.get_wrapper(prefix+fieldName).value()

        return field.to_numpy()


    def getElementCenter(self, filterGhost=False, **kwargs):
        """
        Get element center position as numpy array

        Returns
        -------
            elementCenter : array-like
                Element center coordinates
        """
        elementCenter = self.getField("elementCenter", **kwargs)
        
        if filterGhost:
            elementCenter = self.filterGhostRank(elementCenter, **kwargs)

        return elementCenter

        
    def getElementCenterZ(self, **kwargs):
        """
        Get the z coordinate of the element center

        Returns
        -------
            elementCenterZ : array-like
                Element center z coordinates
        """        
        elementCenter = self.getField("elementCenter", **kwargs)
        elementCenterZ = np.ascontiguousarray(elementCenter[:,2])

        return elementCenterZ


    def getGhostRank(self, **kwargs):
        """
        Get the local ghost ranks
        
        Returns
        -------
            ghostRank : array-like
                Local ghost ranks
        """
        ghostRank = self.getField("ghostRank", **kwargs)

        return ghostRank


    def getLocalToGlobalMap(self, filterGhost=False, **kwargs):
        """
        Get the local rank element id list

        Returns
        -------
            Numpy Array : Array containing the element id list for the local rank
        """
        localToGlobalMap = self.getField("localToGlobalMap", **kwargs)
        
        if filterGhost:
            localToGlobalMap = self.filterGhostRank(localToGlobalMap, **kwargs)
            
        return localToGlobalMap


    def gatherField(self, field, comm, root=0, **kwargs):
        """
        Gather a full GEOS field from all local ranks

        Parameters
        -----------
            field : numpy array
                Local field
            comm : MPI.COMM_WORLD
                MPI communicator
            root : int
                MPI rank used for the gather \
                Default is rank 0
        """
        assert isinstance(root, int)
        assert root < comm.Get_size()

        rank = comm.Get_rank()

        ghostRank = self.getGhostRank(**kwargs)
        localToGlobalMap = self.getLocalToGlobalMap(**kwargs)

        # Prepare buffer 
        nlocalElements = ghostRank.shape[0]
        nmax = np.zeros( 1 )
        nmax[0] = np.max( localToGlobalMap ) # max global number of elements
 
        comm.Barrier()
        comm.Allreduce( MPI.IN_PLACE, nmax, op=MPI.MAX )
        ntot = round( nmax[0] + 1 )

        if rank != root: 
            fullField = None
            nrcv = nlocalElements
            comm.send(nrcv, dest=root, tag=1)
            comm.Send(field, dest=root, tag=2)
            comm.Send(ghostRank, dest=root, tag=3)
            comm.Send(localToGlobalMap, dest=root, tag=4)

        else:
            fullField = np.full( (ntot), fill_value=np.nan )
            jj = np.where( ghostRank < 0 )[0]
            fullField[localToGlobalMap[jj]] = field[jj]

            for r in range( comm.Get_size() ):
                if r != root:
                    nrcv = comm.recv(source=r, tag=1)
    
                    fieldRcv = np.zeros(nrcv, dtype=np.float64)
                    ghostRankRcv = np.zeros(nrcv, dtype=np.int32)
                    localToGlobalMapRcv = np.zeros(nrcv, dtype=np.int64)
    
                    comm.Recv(fieldRcv, source=r, tag=2)
                    comm.Recv(ghostRankRcv, source=r, tag=3)
                    comm.Recv(localToGlobalMapRcv, source=r, tag=4)
        
                    jj = np.where ( ghostRankRcv < 0 )[0]

                    fullField[ localToGlobalMapRcv[jj]] = fieldRcv[jj]                    
        comm.Barrier()
        return fullField, ntot


    def bcastField( self, fullField, comm, root=0, **kwargs ):
        """
        Broadcast a field to local ranks with GEOS local to global map

        Parameters
        -----------
            fullField : numpy array
                Full field
            comm : MPI.COMM_WORLD
                MPI communicator
            root : int
                MPI rank used for the gather \
                Default is rank 0

        Returns
        --------
            field : numpy array
                Local field
        """
        rank = comm.Get_rank()
        size = comm.Get_size()

        ghostRank = self.getGhostRank( **kwargs )
        localToGlobalMap = self.getLocalToGlobalMap( **kwargs )
        nlocalElements = ghostRank.shape[0]

        field = np.zeros( nlocalElements )

        if rank == root:
            jj = np.where(ghostRank < 0)[0]
            field[jj] = fullField[localToGlobalMap[jj]]

            for r in range( size ):
                if r != root:
                    nrcv = comm.recv( source=r, tag=1 )
                    fieldRcv = np.zeros( nrcv, dtype=np.float64 )
                    ghostRankRcv = np.zeros( nrcv, dtype=np.int32 )
                    localToGlobalMapRcv = np.zeros( nrcv, dtype=np.int64 )
    
                    comm.Recv( ghostRankRcv, r, 3 )
                    comm.Recv( localToGlobalMapRcv, r, 4 )
    
                    jj = np.where(ghostRankRcv < 0)[0]
                    fieldRcv[jj] = fullField[localToGlobalMapRcv[jj]]
    
                    comm.Send( fieldRcv, dest=r, tag=100+r )

        else:
            nrcv = nlocalElements
            comm.send( nrcv, root, 1 )
            comm.Send( ghostRank, root, 3 )
            comm.Send( localToGlobalMap, root, 4 )

            comm.Recv( field, source=root, tag=100+rank )

        return field 


    def filterGhostRank(self, field, **kwargs):
        """
        Filter the ghost rank from a GEOS field

        Parameters
        -----------
            field : numpy array
                Field to filter

        Returns
        -------
            field : numpy array
                Filtered field
        """
        ghostRank = self.getGhostRank(**kwargs)
        ind = np.where(ghostRank<0)[0]
        
        return field[ind]


    def getWrapper(self, path):
        """
        Get the GEOS wrapper

        Parameters
        -----------
            path : str
                GEOS path

        Returns
        --------
           Requested wrapper 
        """
        if hasattr(self, "solver"):
            wrapper = self.solver.get_wrapper(path)

            return wrapper


    def getGroup(self, path):
        """
        Get the GEOS group

        Parameters
        -----------
            path : str
                GEOS path

        Returns
        --------
            Group of the path requested
        """
        if hasattr(self, "solver"):
            group = self.solver.get_group(path)

            return group

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
import shutil
import numpy as np
import h5py

import pygeosx
from .WaveSolver import WaveSolver



class AcousticSolver(WaveSolver):
    """
    AcousticSolver Object containing all methods to run AcousticSEM simulation with GEOSX

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


    def setModelForGradient(self, model):
        """
        Set the model for the gradient

        Parameters
        -----------
            model : str
                Model for the velocity
                'c', '1/c' or '1/c2'
        """
        self.model = model


    def updateModel( self, filename, low, high, comm, **kwargs ):
        """
        Update velocity model

        Parameters
        -----------
            filename : str
                .hdf5 file where to get the new model
            low : float
                Min value threshold. All new values < low are set to low
            high : float
                Max value threshold. All new values > high are set to high
            comm : MPI_COMM
               MPI communicators
        """
        root = 0
        rank = comm.Get_rank()
        size = comm.Get_size()

        x = None
        if rank == root:
            with h5py.File( filename, 'r' ) as f:
                x = f["velocity"][:]

            imin = np.where( x < low )[0]
            imax = np.where( x > high )[0]
            x[imin] = low
            x[imax] = high

            if self.model == "1/c2":
                x = np.sqrt(1/x)
            elif self.model == "1/c":
                x = 1/x
            elif self.model == "c":
                pass
            else:
                raise ValueError("Not implemented")

        startModel = self.bcastField( x, comm )

        self.updateVelocityModel( startModel )


    def updateVelocityModel(self, vel, **kwargs):
        """
        Update velocity value in GEOS

        Parameters
        ----------
            vel : array
                Values for velocity field
        """
        super().updateVelocityModel(vel, velocityName="acousticVelocity", **kwargs)


    def setConstantVelocityModel(self, vel, velFieldName="acousticVelocity", **kwargs):
        """
        Update velocity value in GEOS, using a constant value.

        Parameters
        ----------
            vel : float
                Value for velocity field
        """
        prefix = self._getPrefixPath(**kwargs)
        
        velocity = self.solver.get_wrapper(prefix+velFieldName).value()
        velocity.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        velocity.to_numpy()[:] = vel


    def getVelocityModel(self, filterGhost=False, **kwargs):
        """
        Get the velocity values

        Parameters
        -----------
            filterGhost : bool
                Filter the ghost ranks

        Returns
        -------
            velocity : numpy array
                Velocity values
        """
        velocity = super().getVelocityModel(velocityName="acousticVelocity", filterGhost=filterGhost, **kwargs)

        return velocity


    def updateDensityModel(self, density, **kwargs):
        """
        Update density values in GEOS
        
        Parameters
        -----------
            density : array
                New values for the density
        """
        super().updateDensityModel( density=density, densityName="acousticDensity", **kwargs )


    def getGradient(self, filterGhost=False, **kwargs):
        """Get the local rank gradient value

        Returns
        --------
            partialGradient : Numpy Array
                Array containing the element id list for the local rank
        """
        partialGradient = self.getField("partialGradient", **kwargs)

        if filterGhost:
            partialGradient = self.filterGhostRank(partialGradient, **kwargs)

        return partialGradient


    def computePartialGradient( self, shotId, minDepth, comm, gradDirectory="partialGradient", **kwargs ):
        """Compute the partial Gradient

        Parameters
        -----------
            shotId : string
                Number of the shot as string
            minDepth : float
                Depth at which gradient values are kept, otherwise it is set to 0.
                NOTE : this is done in this routine to avoid storage \
                    of elementCenter coordinates in the .hdf5 \
                    but might be problem for WolfeConditions later on \
                    if minDepth is too large
            comm : MPI_COMM
                MPI communicators
            gradDirectory : str, optional
                Partial gradient directory \
                Default is `partialGradient`
        """
        rank = comm.Get_rank()
        size = comm.Get_size()
        root = 0
    
        # Get local gradient
        grad = self.getGradient( **kwargs )
        if self.model == "1/c2":
            x = self.getVelocityModel( **kwargs )
            grad = -( x * x * x / 2 ) * grad
        elif self.model == "1/c":
            x = self.getVelocityModel( filterGhost=True, **kwargs )
            grad = - x * x * grad
        elif self.model == "c":
            pass
        else:
            raise ValueError("Not implemented")

        grad = grad.astype( np.float64 )

        zind = np.where(self.getElementCenter(**kwargs)[:,2] < minDepth)[0]
        grad[zind] = 0.0

        # Gather global gradient
        gradFull, ntot = self.gatherField( field=grad, comm=comm, root=root )        

        if rank == root: 
            os.makedirs(gradDirectory, exist_ok=True)

            with h5py.File(f"{gradDirectory}/partialGradient_"+shotId+".hdf5", 'w') as h5p:
                
                h5p.create_dataset("partialGradient", data=np.zeros(ntot), chunks=True, maxshape=(ntot,))
                h5p["partialGradient"][:] = self.dtWaveField * gradFull
    
            shutil.move(f"{gradDirectory}/partialGradient_"+shotId+".hdf5", f"{gradDirectory}/partialGradient_ready_"+shotId+".hdf5")

        comm.Barrier()


    def getPressureAtReceivers(self):
        """
        Get the pressure values at receivers coordinates

        Returns
        ------
            numpy Array : Array containing the pressure values at all time step at all receivers coordinates
        """
        pressureAtReceivers = self.solver.get_wrapper("pressureNp1AtReceivers").value()
        
        return pressureAtReceivers.to_numpy()
    
    
    def getFullPressureAtReceivers(self, comm):
        """Return all pressures at receivers values on all ranks
        Note that for a too large 2d array this may not work.

        Parameters:
        -----------
            comm : MPI_COMM
                MPI communicators
        """
        rank = comm.Get_rank()

        allPressure = comm.gather(self.getPressureAtReceivers(), root=0)
        pressure = np.zeros(self.getPressureAtReceivers().shape)

        if rank == 0:
            for p in allPressure:
                for i in range(p.shape[1]):
                    if any(p[1:, i])==True:
                        pressure[:, i] = p[:, i]

        pressure = comm.bcast(pressure, root=0)

        return pressure


    def resetWaveField(self, **kwargs):
        """Reinitialize all pressure values on the Wavefield to zero in GEOSX"""

        self.geosx.get_wrapper("Solvers/"+self.name+"/indexSeismoTrace").value()[0] = 0
        meshName = self._getMeshName()
        discretization = self._getDiscretization()

        nodeManagerPath = f"domain/MeshBodies/{meshName}/meshLevels/{discretization}/nodeManager/"


        if self.type == "AcousticSEM":
            for ts in ("nm1", "n", "np1"):
                pressure = self.geosx.get_wrapper(nodeManagerPath + f"pressure_{ts}").value()
                pressure.set_access_level(pygeosx.pylvarray.MODIFIABLE)

                pressure.to_numpy()[:] = 0.0

        elif self.type == "AcousticFirstOrderSEM":
            pressure_np1 = self.geosx.get_wrapper(nodeManagerPath + "pressure_np1").value()
            pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

            pressure_np1.to_numpy()[:] = 0.0

            prefix = self._getPrefixPath(**kwargs)
            for component in ("x", "y", "z"):
                velocity = self.geosx.get_wrapper(prefix + f"velocity_{component}").value()
                velocity.set_access_level(pygeosx.pylvarray.MODIFIABLE)
                
                velocity.to_numpy()[:]   = 0.0


    def resetPressureAtReceivers(self):
        """Reinitialize pressure values at receivers to 0
        """
        pressure = self.solver.get_wrapper("pressureNp1AtReceivers").value()
        pressure.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        pressure.to_numpy()[:] = 0.0


    def getWaveField(self):
        return self.getPressureAtReceivers()[:,:-1]


    def getFullWaveFieldAtReceivers(self, comm):
        return self.getFullPressureAtReceivers(comm)[:,:-1]

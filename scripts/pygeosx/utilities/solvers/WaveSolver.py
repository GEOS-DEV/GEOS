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

from scipy.fftpack import fftfreq, ifft, fft
import numpy as np
import pygeosx

from .Solver import Solver

class WaveSolver(Solver):
    """
    WaveSolver Object containing methods useful for simulation using wave solvers in GEOS

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
        super().__init__(**kwargs)

        self.name = None

        self.dt = dt
        self.minTime = minTime
        self.maxTime = maxTime

        self.minTimeSim = minTime
        self.maxTimeSim = maxTime

        self.dtSeismo = dtSeismo

        self.sourceType = sourceType 
        self.sourceFreq = sourceFreq

        self.dtWaveField = dtWaveField
        

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


    def _setTimeVariables(self):
        """Initialize the time variables"""
        if hasattr(self.xml, "events"):
            events = self.xml.events
            if self.dt is None:
                for event in events["PeriodicEvent"]:
                    if isinstance(event, dict):
                        if event["target"] == "/Solvers/" + self.name:
                            self.dt = float(event['forceDt'])
                    else:
                        if event == "target" and events["PeriodicEvent"]["target"] == "/Solvers/" + self.name:
                            self.dt = float(events["PeriodicEvent"]["forceDt"])
            else:
                self.updateTimeVariable("dt")
 
            if self.maxTime is None:
                self.maxTime = float(events["maxTime"])
                self.maxTimeSim = self.maxTime
            else:
                self.updateTimeVariable("maxTime")
                
            if self.dtSeismo is None:
                if self.type is None:
                    self._getType()
                
                solverdict = self.xml.solvers[self.type]
                for k, v in solverdict.items():
                    if k == "dtSeismoTrace":
                        self.dtSeismo = float(v)
            else:
                self.updateTimeVariable("dtSeismo")

            assert None not in (self.dt, self.dtSeismo, self.maxTime)


    def _setSourceProperties(self):
        """Set the frequency and type of source"""
        if self.sourceFreq is None:
            if self.type is None:
                self._getType()
            solverdict = self.xml.solvers[self.type]
            for k, v in solverdict.items():
                if k == "timeSourceFrequency":
                    self.sourceFreq = v
        
        if self.sourceType is None:        
            if hasattr(self.xml, "events"):
                events = self.xml.events
            try:
                for event in events["PeriodicEvent"]:
                    if isinstance(event, dict):
                        if event["target"] == "/Solvers/" + self.name:
                            self.sourceType = "ricker" + event['rickerOrder']
                    else:
                        if event == "target" and events["PeriodicEvent"]["target"] == "/Solvers/" + self.name:
                            self.sourceType = "ricker" + events["PeriodicEvent"]["rickerOrder"]
            except:
                self.sourceType = "ricker2"


    def initialize(self, rank=0, xml=None):
        """
        Initialization or reinitialization of GEOSX

        Parameters
        ----------
            rank : int
                Process rank
            xml : XML
                XML object containing parameters for GEOSX initialization.
                Only required if not set in the __init__ OR if different from it
        """
        super().initialize(rank, xml)
        self._setSourceProperties()
                

    def updateTimeStep(self, dt):
        """
        Update the solver time step

        Parameters
        ----------
            dt : double
                Time step
        """
        self.dt = dt


    def updateTimeVariable(self, timeVariable):
        """
        Overwrite a GEOS time variable

        Parameters
        ----------
            timeVariable : str
                Name of the time variable to update
        """
        if timeVariable == "maxTime":
            assert hasattr(self, "maxTime")
            self.geosx.get_wrapper("Events/maxTime").value()[0] = self.maxTime
    
        elif timeVariable == "minTime":
            assert hasattr(self, "minTime")
            self.geosx.get_wrapper("Events/minTime").value()[0] = self.minTime

        elif timeVariable == "dt":
            assert hasattr(self, "dt")
            self.geosx.get_wrapper("Events/solverApplications/forceDt").value()[0] = self.dt

        elif timeVariable == "dtSeismo":
            assert hasattr(self, "dtSeismo")
            self.geosx.get_wrapper("/Solvers/"+self.name+"/dtSeismoTrace").value()[0] = self.dtSeismo
 

    def updateSourceFrequency(self, freq):
        """
        Overwrite GEOSX source frequency
        
        Parameters
        ----------
            freq : float
                Frequency of the source in Hz
        """
        self.geosx.get_wrapper("/Solvers/"+self.name+"/timeSourceFrequency").value()[0] = freq
        self.sourceFreq = freq

      

    def updateSourceAndReceivers(self, sourcesCoords=[], receiversCoords=[]):
        """
        Update sources and receivers positions in GEOS

        Parameters
        ----------
            sourcesCoords : list
                List of coordinates for the sources
            receiversCoords : list
                List of coordinates for the receivers
        """
        src_pos_geosx = self.solver.get_wrapper( "sourceCoordinates" ).value()
        src_pos_geosx.set_access_level( pygeosx.pylvarray.RESIZEABLE )
        
        rcv_pos_geosx = self.solver.get_wrapper( "receiverCoordinates" ).value()
        rcv_pos_geosx.set_access_level( pygeosx.pylvarray.RESIZEABLE )

        src_pos_geosx.resize_all(( len( sourcesCoords ), 3 ))
        if len( sourcesCoords ) == 0:
            src_pos_geosx.to_numpy()[:] = np.zeros(( 0, 3 ))
        else:
            src_pos_geosx.to_numpy()[:] = sourcesCoords[:]


        rcv_pos_geosx.resize_all(( len( receiversCoords ), 3 ))
        if len( receiversCoords ) == 0:
            rcv_pos_geosx.to_numpy()[:] = np.zeros(( 0, 3 ))
        else:
            rcv_pos_geosx.to_numpy()[:] = receiversCoords[:]

        self.solver.reinit()


    def evaluateSource(self):
        """
        Evaluate source and update on GEOS
        Only ricker order {0 - 4} accepted
        """
        sourceTypes = ("ricker0", "ricker1", "ricker2", "ricker3", "ricker4")
        assert self.sourceType in sourceTypes, f"Only {sourceTypes} are allowed"

        f0 = self.sourceFreq
        delay = 1.0 / f0
        alpha = - ( f0 * np.pi )**2

        nsamples = int( round( ( self.maxTime - self.minTime ) / self.dt )) + 1
        sourceValue = np.zeros(( nsamples, 1 ))
        
        order = int( self.sourceType[-1] )
        sgn = ( -1 )**( order + 1 )

        time = self.minTime
        for nt in range(nsamples):

            if self.minTime <= - 1.0 / f0:
                tmin = - 2.9 / f0
                tmax = 2.9 / f0
                time_d = time

            else:
                time_d = time - delay
                tmin = 0.0
                tmax = 2.9 / f0

            if (time > tmin and time < tmax) or ( self.minTime < - 1 / f0 and time == tmin ):
                gaussian = np.exp( alpha * time_d**2)

                if order == 0:
                    sourceValue[nt, 0] = sgn * gaussian

                elif order == 1:
                    sourceValue[nt, 0] = sgn * ( 2 * alpha * time_d ) * gaussian

                elif order == 2:
                    sourceValue[nt, 0] = sgn * ( 2 * alpha + 4 * alpha**2 * time_d**2 ) * gaussian

                elif order == 3:
                    sourceValue[nt, 0] =  sgn * ( 12 * alpha**2 * time_d + 8 * alpha**3 * time_d**3 ) * gaussian

                elif order == 4:
                    sourceValue[nt, 0] = sgn * ( 12 * alpha**2  + 48 * alpha**3 * time_d**2 + 16 * alpha**4 * time_d**4 ) * gaussian

            time += self.dt
         
        self.updateSourceFrequency(self.sourceFreq)
        self.updateSourceValue(sourceValue)
        self.sourceValue = sourceValue


    def updateSourceValue(self, value):
        """
        Update the value of the source in GEOS

        Parameters
        ----------
            value : array/list
                List/array containing the value of the source at each time step
        """
        src_value = self.solver.get_wrapper("sourceValue").value()
        src_value.set_access_level(pygeosx.pylvarray.RESIZEABLE)

        src_value.resize_all(value.shape)
        src_value.to_numpy()[:] = value[:]

        self.maxTimeSim = ( value.shape[0] - 1 ) * self.dt
        self.geosx.get_wrapper("Events/minTime").value()[0] = self.minTime
        self.sourceValue = value[:]


    def filterSource(self, fmax):
        """
        Filter the source value and give the value to GEOSX. Note that is can also modify the start and end time of simulation to avoid discontinuity.

        Parameters
        -----------
            fmax : float/string
                Max frequency of the source wanted. The source then have frequencies in the interval [0,fmax+1]
        """
        if str(fmax) == "all":
            return

        minTime = self.minTime
        maxTime = self.maxTime
        dt = self.dt
        f0 = self.sourceFreq

        sourceValue = self.sourceValue

        pad = int(round(sourceValue.shape[0]/2))
        n = sourceValue.shape[0] + 2 * pad

        tf = fftfreq(n, dt)
        y_fft = np.zeros((n,sourceValue.shape[1]), dtype="complex_")
        y = np.zeros(y_fft.shape, dtype="complex_")

        for i in range(y_fft.shape[1]):
            y_fft[pad:n-pad,i] = sourceValue[:,i]
            y_fft[:,i] = fft(y_fft[:,i])# Perform fourier transform

        isup = np.where(tf>=fmax)[0]
        imax = np.where(tf[isup]>=fmax+1)[0][0]
        i1 = isup[0]
        i2 = isup[imax]

        iinf = np.where(tf<=-fmax)[0]
        imin = np.where(tf[iinf]<=-fmax-1)[0][-1]

        i3 = iinf[imin]
        i4 = iinf[-1]

        for i in range(y_fft.shape[1]):
            y_fft[i1:i2,i] = np.cos((isup[0:imax] - i1)/(i2-i1) * np.pi/2)**2 * y_fft[i1:i2,i]
            y_fft[i3:i4,i] = np.cos((iinf[imin:-1] - i4)/(i3-i4) * np.pi/2)**2 * y_fft[i3:i4,i]
            y_fft[i2:i3,i] = 0

        for i in range(y.shape[1]):
            y[:,i] = ifft(y_fft[:,i])# Perform inverse fourier transform

        it0 = int(round(abs(minTime/dt))) + pad
        d = int(round(1/f0/dt))

        i1 = max(it0 - 4*d, 0)
        i2 = int(round(i1 + d/4))

        i4 = min(n,n - pad + 4*d)
        i3 = int(round(i4 - d/4))

        for i in range(y.shape[1]):
            y[i1:i2,i] = np.cos((np.arange(i1,i2) - i2)/(i2-i1) * np.pi/2)**2 * y[i1:i2,i]
            y[i3:i4,i] = np.cos((np.arange(i3,i4) - i3)/(i4-i3) * np.pi/2)**2 * y[i3:i4,i]
            y[max(i1-d,0):i1,i] = 0.0
            y[i4:min(i4+d,n),i] = 0.0


        t = np.arange(minTime-pad*dt, maxTime+pad*dt+dt/2, dt)

        self.updateSourceValue(np.real(y[max(i1-d,0):min(i4+d,n),:]))
        self.minTimeSim = t[max(i1-d,0)]
        self.maxTimeSim = t[min(i4+d,n-1)]
        self.geosx.get_wrapper("Events/minTime").value()[0] = self.minTimeSim
        self.sourceValue = np.real(y[max(i1-d,0):min(i4+d,n),:])


    def updateVelocityModel(self, vel, velocityName, **kwargs):
        """
        Update velocity value in GEOS

        Parameters
        ----------
            vel : float/array
                Value(s) for velocity field
            velocityName : str
                Name of the velocity array in GEOS
        """
        prefix = self._getPrefixPath(**kwargs)
        
        velocity = self.solver.get_wrapper(prefix+velocityName).value()
        velocity.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        velocity.to_numpy()[vel>0] = vel[vel>0]


    def updateDensityModel( self, density, densityName, **kwargs ):
        """
        Update density values in GEOS
        
        Parameters
        -----------
            density : array
                New values for the density
            densityName : str
                Name of density array in GEOS
        """
        prefix = self._getPrefixPath( **kwargs )

        d = self.solver.get_wrapper( prefix + densityName ).value()
        d.set_access_level( pygeosx.pylvarray.MODIFIABLE )

        d.to_numpy()[density>0] = density[density>0]


    def outputWaveField(self, time):
        """
        Trigger the wavefield output

        Parameters
        ----------
            time : float
                Current time of simulation
        """
        self.collections[0].collect(time, self.dt)
        self.hdf5Outputs[0].output(time, self.dt)


    def getVelocityModel(self, velocityName, filterGhost=False, **kwargs):
        """
        Get the velocity values

        Parameters
        -----------
            velocityName : str
                Name of velocity array in GEOS
            filterGhost : bool
                Filter the ghost ranks

        Returns
        -------
            Numpy Array : Array containing the velocity values
        """
        velocity = self.getField(velocityName, **kwargs)
        
        if filterGhost:
            velocity = self.filterGhostRank(velocity, **kwargs)

        return velocity


    def getWaveField(self):
        pass

    def getWaveFieldAtReceivers(self, comm):
        pass

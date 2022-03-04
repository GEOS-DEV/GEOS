import pygeosx
import numpy as np
import os
import xml.etree.ElementTree as ET
import sys

class AcousticSolver:
    def __init__(self,
                 xml,
                 dt=None,
                 maxTime=None,
                 dtSeismo=None,
                 **kwargs):

        """Creates AcousticSolver Object that contains all methods to run AcousticSEM simulation with GEOSX

        Parameters
        ----------
        xml : xml file
            XML file containing parameters for GEOSX initialization

        dt : float
            Time step for simulation

        maxTime : float
            End Time of simulation

        dtSeismo : float
            Time step to save pressure for seismic trace
        """

        self.xml = xml

        tree = ET.parse(self.xml)
        root = tree.getroot()

        acousticSEM = [elem.attrib for elem in root.iter('AcousticSEM')]
        events = [elem.attrib for elem in root.iter('PeriodicEvent')]
        hdf5s = [elem.attrib for elem in root.iter('TimeHistory')]

        if len(acousticSEM):
            self.name = acousticSEM[0]['name']
        else:
            raise ValueError("You must set AcousticSEM section in XML.")

        if dt is not None:
            self.dt = dt
        else:
            for event in events:
                if event['target'] == "/Solvers/"+self.name:
                    self.dt = float(event['forceDt'])

        if maxTime is not None:
            self.maxTime = maxTime
        else:
            maxTime = [elem.attrib for elem in root.iter('Events')][0]['maxTime']
            self.maxTime = maxTime

        if dtSeismo is not None:
            self.dtSeismo = dtSeismo
        else:
            for key in acousticSEM.keys():
                if key == "dtSeismoTrace":
                    self.dtSeismo = acousticSEM[key]

        self.collectionTargets = []
        self.outputTargets = []
        for hdf5 in hdf5s:
            self.collectionTargets.append(hdf5['sources'][1:-1])
            self.outputTargets.append("Outputs/"+hdf5['name'])


    def __repr__(self):
        string_list = []
        string_list.append("Solver type : " + type(self).__name__ + "\n")
        string_list.append("dt : " + str(self.dt) + "\n")
        string_list.append("maxTime : " + str(self.maxTime) + "\n")
        string_list.append("dtSeismo : " + str(self.dtSeismo) +"\n")
        string_list.append("Outputs : " + str(self.outputTargets) +"\n")
        rep=""
        for string in string_list:
            rep += string

        return rep


    def initialize(self, rank=0):

        """First initialization of GEOSX.
        If you have already run GEOSX, use the reinit() method instead

        Parameters
        ----------
        rank : int
            Process rank
        """

        self.geosx = pygeosx.initialize(rank, sys.argv)

        self.solver = self.geosx.get_group("/Solvers/"+self.name)
        self.collections = []
        self.outputs = []

        for target in self.collectionTargets:
            self.collections.append(self.geosx.get_group(target))

        for target in self.outputTargets:
            self.outputs.append(self.geosx.get_group(target))

        self.updateTimeVariables()


    def reinit(self, rank=0):

        """First initialization of GEOSX.
        Must be used only if you have already run GEOSX once

        Parameters
        ----------
        rank : int
            Process rank
        """

        self.geosx = pygeosx.reinit(rank)

        self.solver = self.geosx.get_group("/Solvers/"+self.name)
        self.collections = []
        self.outputs = []

        for target in self.collectionTargets:
            self.collections.append(self.geosx.get_group(target))

        for target in self.outputTargets:
            self.outputs.append(self.geosx.get_group(target))

        self.updateTimeVariables()


    def apply_initial_conditions(self):

        """Second initialization of GEOSX.
        """

        pygeosx.apply_initial_conditions()


    def finalize(self):

        """Terminate GEOSX
        """

        pygeosx._finalize()


    def updateTimeVariables(self):

        """Overwrite GEOSX time variables that have been read in the XML.
        """

        self.geosx.get_wrapper("Events/maxTime").value()[0] = self.maxTime
        self.geosx.get_wrapper("Events/solverApplications/forceDt").value()[0] = self.dt
        self.geosx.get_wrapper("/Solvers/"+self.name+"/dtSeismoTrace").value()[0] = self.dtSeismo


    def updateOutputsName(self, list_of_output):

        """Overwrite GEOSX hdf5 Outputs paths that have been read in the XML.

        Parameters
        ----------
        list_of_output : list
            List of string path
        """

        if not len(self.outputs):
            raise ValueError("No Outputs specified in XML.")
        else:
            for i in range(len(list_of_output)):
                if os.path.exists(os.path.dirname(list_of_output[i])):
                    pass
                else:
                    os.system("mkdir -p " + list_of_output[i])

                self.outputs[i].setOutputName(list_of_output[i])
                self.outputs[i].reinit()


    def updateSourceAndReceivers(self, sources_list=[], receivers_list=[]):

        """Update Sources and Receivers positions in GEOSX.

        Parameters
        ----------
        sources_list : list
            List of coordinates for the sources

        receivers_list : list
            List of coordinates for the receivers
        """

        src_pos_geosx = self.solver.get_wrapper("sourceCoordinates").value()
        src_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)

        rcv_pos_geosx = self.solver.get_wrapper("receiverCoordinates").value()
        rcv_pos_geosx.set_access_level(pygeosx.pylvarray.RESIZEABLE)


        src_pos_geosx.resize(len(sources_list))
        if len(sources_list) == 0:
            src_pos_geosx.to_numpy()[:] = np.zeros((0,3))
        else:
            src_pos = [source.coords for source in sources_list]
            src_pos_geosx.to_numpy()[:] = src_pos[:]

        rcv_pos_geosx.resize(len(receivers_list))
        if len(receivers_list) == 0:
            rcv_pos_geosx.to_numpy()[:] = np.zeros((0,3))
        else:
            rcv_pos = [receiver.coords for receiver in receivers_list]
            rcv_pos_geosx.to_numpy()[:] = rcv_pos[:]

        self.solver.reinit()


    def updateSourceValue(self, value):

        """Update the value of the source in GEOSX

        Parameters
        ----------
        value : array/list
            List/array containing the value of the source at each time step
        """

        src_value = self.solver.get_wrapper("sourceValue").value()
        src_value.set_access_level(pygeosx.pylvarray.MODIFIABLE)
        src_value.to_numpy()[:] = value[:]


    def resetWaveField(self):

        """Reinitialize all pressure values on the Wavefield to 0.
        """

        self.geosx.get_wrapper("Solvers/"+self.name+"/indexSeismoTrace").value()[0] = 0
        nodeManagerPath = "domain/MeshBodies/mesh/meshLevels/Level0/nodeManager/"

        pressure_nm1 = self.geosx.get_wrapper(nodeManagerPath + "pressure_nm1").value()
        pressure_nm1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        pressure_n = self.geosx.get_wrapper(nodeManagerPath + "pressure_n").value()
        pressure_n.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        pressure_np1 = self.geosx.get_wrapper(nodeManagerPath + "pressure_np1").value()
        pressure_np1.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        pressure_geosx = self.geosx.get_wrapper("Solvers/"+self.name+"/pressureNp1AtReceivers").value()
        pressure_geosx.set_access_level(pygeosx.pylvarray.MODIFIABLE)

        pressure_nm1.to_numpy()[:] = 0.0
        pressure_n.to_numpy()[:]   = 0.0
        pressure_np1.to_numpy()[:] = 0.0
        pressure_geosx.to_numpy()[:] = 0.0


    def execute(self, time):

        """Do one iteration of the solver.

        Parameters
        ----------
        time : float
            Current time of simulation
        """

        self.solver.execute(time, self.dt)


    def outputWaveField(self, time):

        """Trigger the output

        Parameters
        ----------
        time : float
            Current time of simulation
        """

        for i in range(len(self.outputs)):
            self.collections[i].collect(time, self.dt)
            self.outputs[i].output(time, self.dt)


    def getPressureAtReceivers(self):

        """Get the pressures values at receivers coordinates

        Return
        ------
        pressureAtReceivers : Numpy Array
            Array containing the pressure values at all time step at all receivers coordinates
        """

        pressureAtReceivers = self.solver.get_wrapper("pressureNp1AtReceivers").value()
        return pressureAtReceivers.to_numpy()

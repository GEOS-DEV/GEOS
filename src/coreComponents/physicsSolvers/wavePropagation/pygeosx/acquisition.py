from copy import deepcopy
import numpy as np
import segyio
import random
import os


class Acquisition:
    def __init__(self,
                 limited_aperture=False,
                 aperture_boundary=None,
                 dt=None,
                 velocity_model=None,
                 acquisition_type=None,
                 boundary=None,
                 source=None,
                 source_positions=None,
                 receivers_positions=None,
                 source_depth=None,
                 receivers_depth=None,
                 output=None):

        self.limited_aperture=limited_aperture
        self.boundary=boundary

        if source_positions is not None and receivers_positions is not None:
            self.acquisition(source,
                             source_positions,
                             receivers_positions,
                             source_depth,
                             receivers_depth)


        dist_from_source = 500
        depth_from_source = 500
        if limited_aperture is False:
            self.velocity_model=velocity_model
            self.aperture_boundary=boundary
            if dt is not None:
                self.dt=dt
            else:
                self.dt=calculDt(velocity_model,
                                 boundary)
        else:
            self.velocity_model=[]

            self.set_limited_aperture_boundaries(dist_from_source, depth_from_source)
            for i in range(len(self.shots)):
                self.shots[i].receivers.keep_inside_domain(self.limited_aperture[i])


        if output is not None:
            self.output=output
        else:
            self.output = os.path.join(os.getcwd(),"outputSismoTrace/")
            if os.path.exists(self.output):
                pass
            else:
                os.mkdir(self.output)


    def acquisition(self,
                    wavelet,
                    source_positions,
                    receivers_positions,
                    source_depth,
                    receivers_depth):

        nr=len(receivers_positions)
        ns=len(source_positions)

        receivers = ReceiverSet([Receiver([receivers_positions[i,0], receivers_positions[i,1], receivers_depth])
                                 for i in range(nr)])

        shots = []

        for i in range(ns):

            if len(str(i+1))<2:
                shot_id = "00"+str(i+1)
            elif len(str(i+1))<3 and len(str(i+1))>2:
                shot_id = "0"+str(i+1)
            else:
                shot_id = str(i+1)

            # Define source location and type
            source = Source([source_positions[i,0], source_positions[i,1], source_depth], wavelet)

            shot = Shot(source, receivers, shot_id)
            shots.append(shot)

        self.shots = shots



    def split(self, n):
        acqs = []
        nb_shot_m1 = len(self.shot)
        ind = 0

        for i in range(n):
            nb_shot=int(nb_shot_m1/(n-i))

            acqs.append(deepcopy(self))
            acqs[i].shot = acqs[i].shot[ind:ind + nb_shot]

            ind = ind + nb_shot
            nb_shot_m1 = nb_shot_m1 - nb_shot

        return acqs


    def setDt(self, dt):
        self.dt=dt



    def set_limited_aperture_boundaries(self, distance, depth):
        self.aperture_boundary=[]

        for shot in self.shots:
            if shot.source.x() - distance < self.boundary[0][0]:
                xmin = boundary[0][0]
            else:
                xmin = shot.source.x() - distance
            if shot.source.x() + distance > self.boundary[0][1]:
                xmax = boundary[0][1]
            else:
                xmax = shot.source.x() + distance

            if shot.source.y() - distance < self.boundary[1][0]:
                ymin = boundary[1][0]
            else:
                ymin = shot.source.y() - distance
            if shot.source.y() + distance > self.boundary[1][1]:
                ymax = boundary[1][1]
            else:
                ymax = shot.source.y() + distance

            zmin = boundary[2][0]

            if shot.source.z() + depth > self.boundary[2][1]:
                zmax = boundary[2][1]
            else:
                zmax = shot.source.z() + depth

            limited_aperture = [[xmin,xmax],[ymin,ymax],[zmin,zmax]]
            self.aperture_boundary.append(copy.deepcopy(limited_aperture))



    def construct_from_dict(self, **kwargs):
        for key, value in kwargs.items():
            if key == "shot":
                setattr(self, key, Shot().construct_from_dict(**value))
            else:
                setattr(self, key, value)

        return self




class SEGYAcquisition(Acquisition):

    def __init__(self,
                 limited_aperture=False,
                 velocity_model=None,
                 boundary=None,
                 dt=None,
                 directory=None,
                 source=None):

        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model)

        if directory is not None:
            self.acquisition(directory,
                             source)
        else:
            raise ValueError("You must specify a directory")

    def acquisition(self,
                    directory,
                    wavelet):

        """Acquisition from .sgy files

        Parameters
        ----------
        segyPath : string
            Directory with segy files to read

        wavelet : list
            Source function (Ricker)


        Return
        ------
        shots : list of Shot object
        List of shots configuration
        """

        shot_list = []

        i=0
        for filename in glob.glob(os.path.join(directory, '*.sgy')):
            receiver_list = []

            j=0
            with segyio.open(filename, 'r', ignore_geometry=True) as f:

                scalarXY = float(f.header[0][71])
                scalarZ = float(f.header[0][69])

                sourceX = f.header[0][73]*abs(scalarXY)**np.sign(scalarXY)
                sourceY = f.header[0][77]*abs(scalarXY)**np.sign(scalarXY)
                sourceZ = f.header[0][49]*abs(scalarXY)**np.sign(scalarZ)

                source = Source([sourceX, sourceY, sourceZ], wavelet)
                for i in range(len(f.trace)):

                    receiverX = f.header[i][81]*abs(scalarXY)**np.sign(scalarXY)
                    receiverY = f.header[i][85]*abs(scalarXY)**np.sign(scalarXY)
                    receiverZ = f.header[i][41]*abs(scalarXY)**np.sign(scalarZ)

                    receiver_list.append(Receiver([receiverX, receiverY, receiverZ]))

                I = (i+1)*j
                if len(str(I+1))<2:
                    shot_id = "00"+str(I+1)
                elif len(str(I+1))<3 and len(str(I+1))>2:
                    shot_id = "0"+str(I+1)
                else:
                    shot_id = str(I+1)

                receivers = ReceiverSet(receiver_list)

                shot_list.append(Shot(source, receivers, shot_id))
                j+=1
            i+=1

        self.shots = shot_list



class EQUISPACEDAcquisition(Acquisition):

    def __init__(self,
                 limited_aperture=False,
                 velocity_model=None,
                 boundary=None,
                 dt=None,
                 source=None,
                 start_source_pos=None,
                 end_source_pos=None,
                 start_receivers_pos=None,
                 end_receivers_pos=None,
                 number_of_sources=None,
                 number_of_receivers=None,
                 source_depth=None,
                 receivers_depth=None):

        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model,
                         dt=dt)

        if start_source_pos is not None and start_receivers_pos is not None:
            self.acquisition(wavelet=source,
                             start_source_pos=start_source_pos,
                             end_source_pos=end_source_pos,
                             start_receivers_pos=start_receivers_pos,
                             end_receivers_pos=end_receivers_pos,
                             number_of_sources=number_of_sources,
                             number_of_receivers=number_of_receivers,
                             source_depth=source_depth,
                             receivers_depth=receivers_depth)


    def acquisition(self,
                    wavelet=None,
                    start_source_pos=None,
                    end_source_pos=None,
                    start_receivers_pos=None,
                    end_receivers_pos=None,
                    number_of_sources=None,
                    number_of_receivers=None,
                    source_depth=None,
                    receivers_depth=None):

        """This acquisition is for linears and crosses acquisitions (possibility to define
        multiple lines of receivers)

        Parameters
        ----------
        boundary : list
            List containing min/max boundary coordinates of the domain [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

        wavelet : list
            Source function (Ricker)

        start_source_pos : list
            Starting coordinate for source line

        end_source_pos : list
            Ending coordinate for source line

        start_receivers_pos : list of list
            Starting coordinate for receivers line

        end_receivers_pos : list of list
            Ending coordinate for receivers line

        number_of_sources : int
            Number of sources

        number_of_receivers : list
            Number of receivers for the corresponding line (if only one given, number applied to all lines of receivers)

        source_depth : float
            z source position

        receivers_depth : float
            z receivers position


        Return
        ------
        shots : list of Shot object
        list of shots configurations
        """

        receivers = ReceiverSet()

        for i in range(len(start_receivers_pos)):
            if isinstance(number_of_receivers, int):
                xr = np.linspace(start_receivers_pos[i][0], end_receivers_pos[i][0], number_of_receivers).tolist()
                yr = np.linspace(start_receivers_pos[i][1], end_receivers_pos[i][1], number_of_receivers).tolist()
            else:
                xr = np.linspace(start_receivers_pos[i][0], end_receivers_pos[i][0], number_of_receivers[i]).tolist()
                yr = np.linspace(start_receivers_pos[i][1], end_receivers_pos[i][1], number_of_receivers[i]).tolist()

        receivers_temp = ReceiverSet([Receiver([x, y, receivers_depth]) for x, y in list(zip(xr, yr))])
        #print(receivers_temp)
        receivers.append(deepcopy(receivers_temp))

        xs = np.linspace(start_source_pos[0], end_source_pos[0], number_of_sources).tolist()
        ys = np.linspace(start_source_pos[1], end_source_pos[1], number_of_sources).tolist()

        # receivers = receivers.getInsideDomain(box)
        shots = []

        for i in range(len(xs)):

            if len(str(i+1))<2:
                shot_id = "00"+str(i+1)
            elif len(str(i+1))<3 and len(str(i+1))>2:
                shot_id = "0"+str(i+1)
            else:
                shot_id = str(i+1)

            srcpos = [xs[i], ys[i], source_depth]
            source = Source(srcpos, wavelet)

            shot = Shot(source, receivers, shot_id)
            shots.append(shot)

        self.shots = shots


class MOVINGAcquisition(Acquisition):

    def __init__(self,
                 boundary=None,
                 limited_aperture=False,
                 velocity_model=None,
                 dt=None,
                 source=None,
                 number_of_sources_x = None,
                 number_of_sources_y = None,
                 number_of_receivers_x = None,
                 number_of_receivers_y = None,
                 source_depth = None,
                 receivers_depth = None,
                 zone_receiver_x = None,
                 zone_receiver_y = None):

        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model,
                         dt=dt)

        if number_of_sources_x is not None and number_of_receivers_x is not None:
            self.acquisition(boundary=boundary,
                             wavelet=source,
                             start_source_pos=start_source_pos,
                             end_source_pos=end_source_pos,
                             start_receivers_pos=start_receivers_pos,
                             end_receivers_pos=end_receivers_pos,
                             number_of_sources=number_of_sources,
                             number_of_receivers=number_of_receivers,
                             source_depth=source_depth,
                             receivers_depth=receivers_depth)


    def acquisition(self,
                    boundary=None,
                    wavelet=None,
                    number_of_sources_x = None,
                    number_of_sources_y = None,
                    number_of_receivers_x = None,
                    number_of_receivers_y = None,
                    source_depth = None,
                    receivers_depth = None,
                    zone_receiver_x = None,
                    zone_receiver_y = None):

        """Marine seismic acquisition reprensenting a boat pulling a source
        and a set of receivers in a squared domain

        Parameters
        ----------
        boundary : list
            List containing min/max boundary coordinates of the domain [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

        wavelet : list
            Source function (Ricker)

        number_of_sources_x : int
            Number of sources along x axis (this is the number of configuration, 1 position for each configuration)

        number_of_sources_y : int
            Number of sources along y axis

        number_of_receivers_x : int
            Number of sources along x axis

        number_of_receivers_y : int
            Number of sources along y axis

        source_depth : float
            z source position

        receivers_depth : float
            z receivers position

        zone_receiver_x : float
            Lenght of the area for receivers along x axis (distance from the receivers coordinate)

        zone_receiver_y : float
            Lenght of the area for receivers along y axis


        Return
        ------
        shots : list of Shot object
        List of shots configuration
        """


        xmin = boundary[0][0]
        xmax = boundary[0][1]

        ymin = boundary[1][0]
        ymax = boundary[1][1]

        zmin = boundary[2][0]
        zmax = boundary[2][1]

        movex = (xmax-xmin)/(number_of_sources_x + 1) + 1
        movey = (ymax-ymin)/(number_of_sources_y + 1) + 1

        srcpos = [xmin, ymin, source_depth]

        xposleft  = np.linspace(xmin + 0.3*movex, xmin + 0.3*movex + zone_receiver_x, number_of_receivers_x)
        xposright = np.linspace(xmax - 0.3*movex, xmax - 0.3*movex - zone_receiver_x, number_of_receivers_x)

        ypos = np.linspace(ymin - zone_receiver_y/2, ymin + zone_receiver_y/2, number_of_receivers_y).tolist()

        receiversbaseleft  = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xposleft for y in ypos])
        receiversbaseright = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xposright for y in ypos])

        shots = []

        for i in range(number_of_sources_y):
            if i%2==0:
                srcpos[0] = xmin
                receivers = copy.deepcopy(receiversbaseleft)

            else:
                srcpos[0] = xmax
                receivers = copy.deepcopy(receiversbaseright)

            srcpos[1] += movey
            receivers.linearSetMove(1, (i+1)*movey)
            for j in range(number_of_sources_x):

                I = (i+1)*j
                if len(str(I+1))<2:
                    shot_id = "00"+str(I+1)
                elif len(str(I+1))<3 and len(str(I+1))>2:
                    shot_id = "0"+str(I+1)
                else:
                    shot_id = str(I+1)

                if i%2==0:
                    srcpos[0] += movex
                    receivers.linearSetMove(0, movex)

                else:
                    srcpos[0] -= movex
                    receivers.linearSetMove(0, -movex)

                rcvset = receivers.getInsideDomain(boundaries)

                # Define source location and type
                source = Source(srcpos, wavelet)

                shot = Shot(source, copy.deepcopy(rcvset), shot_id)
                shots.append(shot)


        self.shots = shots


class RANDOMAcquisition(Acquisition):

    def __iter__(self,
                 limited_aperture=False,
                 boundary=None,
                 velocity_model=None,
                 dt=None,
                 source=None,
                 number_of_sources = None,
                 source_depth = None,
                 number_of_receivers_x = None,
                 number_of_receivers_y = None,
                 receivers_depth = None):

        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         velocity_model=velocity_model,
                         dt=dt,
                         source=source)

        if number_of_sources is not None and number_of_receivers_x is not None and number_of_receivers_y is not None:
            self.acquisition(wavelet=source,
                             boundary=boundary,
                             number_of_sources=number_of_sources,
                             number_of_receivers_x=number_of_receivers_x,
                             number_of_receivers_y=number_of_receivers_y,
                             source_depth=source_depth,
                             receivers_depth=receivers_depth)


    def acquisition(self,
                    boundary=None,
                    wavelet=None,
                    number_of_sources = None,
                    source_depth = None,
                    number_of_receivers_x = None,
                    number_of_receivers_y = None,
                    receivers_depth = None):

        """Random seismic acquisition, the positions of sources are set randomly,
        the receivers are set as a grid over the domain

        Parameters
        ----------
        boundary : list
            List containing min/max boundary coordinates of the domain [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

        wavelet : list
            Source function (Ricker)

        number_of_sources : int
            Number of sources

        source_depth : float
            z source position

        number_of_receivers_x : int
            Number of sources along x axis

        number_of_receivers_y : int
            Number of sources along y axis

        receivers_depth : float
            z receivers position

        Return
        ------
        shots : list of Shot object
        List of shots configurations
        """

        xmin = boundary[0][0]
        xmax = boundary[0][1]

        ymin = boundary[1][0]
        ymax = boundary[1][1]

        zmin = boundary[2][0]
        zmax = boundary[2][1]


        xpos = np.linspace(xmin, xmax, number_of_receivers_x).tolist()
        ypos = np.linspace(ymin, ymax, number_of_receivers_y).tolist()

        receivers = ReceiverSet([Receiver([x, y, receivers_depth]) for x in xpos for y in ypos])

        shots = list()

        for i in range(number_of_sources):

            if len(str(i+1))<2:
                shot_id = "00"+str(i+1)
            elif len(str(i+1))<3 and len(str(i+1))>2:
                shot_id = "0"+str(i+1)
            else:
                shot_id = str(i+1)

            srcpos = [random.random()*(xmin-0.2*(xmax-xmin)),(xmax+0.2*(xmax-xmin))*random.random(), source_depth]

            # Define source location and type
            source = Source(srcpos, wavelet)

            shot = Shot(source, receivers, shot_id)
            shots.append(shot)

        self.shots = shots




'''Create Shot object composed of one Source and a ReceiverSet'''
class Shot:
    """Class representing a shot configuration

    Attributes
    ----------
    source :
        Source object

    receivers :
        ReceiverSet object

    flag :
        A flag to say if the shot configuration has been simulated
        "Undone", "In Progress", "Done"
    """

    def __init__(self,
                 source=None,
                 receivers=None,
                 shot_id=None):
        """ Constructor of Shot

        Parameters
        ----------
        Source :
            A Source object

        Receivers :
            A ReceiverSet object
        """

        self.source = source
        self.receivers = receivers
        self.flag = "Undone"

        if shot_id is not None:
            self.id = shot_id
        else:
            self.id = None

    def __repr__(self):
        return 'Source position : \n'+str(self.source) +' \n\n' + 'Receivers positions : \n' + str(self.receivers) + '\n\n'


    def flagUpdate(self, string):
        self.flag = string


     def construct_from_dict(self, **kwargs):
        for key, value in kwargs.items():
            if key == "source":
                setattr(self, key, Source().construct_from_dict(**value))
            elif key == "receivers":
                setattr(self, key, ReceiverSet().construct_from_dict(value))
            else:
                setattr(self, key, value)

        return self



class Receiver:
    """A class representing a receiver

    Attributes
    ----------
    coords :
        Coordinates of the source
    """

    def __init__(self, coords=None):
        """Constructor for the receiver

        Parameters
        ----------
        pos :
            Coordinates for the receiver
        """

        self.coords = coords

    def __repr__(self):
        return '('+str(self.coords[0])+','+str(self.coords[1])+','+str(self.coords[2])+') \n'


    def setReceiverPos(self, coords, x, y, z):
        if coord==0:
            self.coords[0] = x
        elif coord==1:
            self.coords[1] = y
        elif coord==2:
            self.coords[2] = z

    #Transpose the receiver position
    def linearMove(self, coord, value):
        if coord==0:
            self.coords[0] += value
        elif coord==1:
            self.coords[1] += value
        elif coord==2:
            self.coords[2] += value

    def x(self):
        return self.coords[0]

    def y(self):
        return self.coords[1]

    def z(self):
        return self.coords[2]


    def construct_from_dict(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return self


class ReceiverSet:
    """A class representing a set receiver

    Attributes
    ----------
    receiver_list :
        List of Receiver

    n :
        Number of Receiver
    """

    def __init__(self, receiver_list = None):
        """Constructor for the point source

        Parameters
        ----------
        receiver_list :
            List of Receiver
        """
        if receiver_list is None:
            self.receiver_list = []
            self.n = 0
        else:
            self.receiver_list = receiver_list
            self.n = len(receiver_list)

    def __repr__(self):
        if self.n >=10:
            return str(self.receiver_list[0:4])[:-1] + '...' + '\n' + str(self.receiver_list[-4:])[1:]
        else:
            return str(self.receiver_list)


    def keep_inside_domain(self, meshBoundaries):
        new_list=[]

        for receiver in self.receiver_list:
            if meshBoundaries[0][0] <= receiver.coords[0] and receiver.coords[0] <= meshBoundaries[0][1] \
            and meshBoundaries[1][0] <= receiver.coords[1] and receiver.coords[1] <= meshBoundaries[1][1] \
            and meshBoundaries[2][0] <= receiver.coords[2] and receiver.coords[2] <= meshBoundaries[2][1]:
                new_list.append(receiver)

        self.receiver_list = new_list
        self.n = len(new_list)


    def getSetCoord(self):
        listRcv = []
        for i in range(self.n):
            listRcv.append(np.array([self.receiver_list[i].coords[0], self.receiver_list[i].coords[1], self.receiver_list[i].coords[2]]))
        return listRcv

    def append(self, ReceiverSet):
        for i in range(ReceiverSet.n):
            self.receiver_list.append(ReceiverSet.receiver_list[i])
        self.n += ReceiverSet.n


    #Transpose all receivers position
    def linearSetMove(self, coord, value):
        for i in range(self.n):
            self.receiver_list[i].linearMove(coord, value)


     def construct_from_dict(self, list_of_dict):
        for value in list_of_dict:
            if key == "receiver":
                setattr(self, key, Receiver().construct_from_dict(**value))
            else:
                setattr(self, key, value)

        return self



class Source:
    """A class representing a point source

    Attributes

    ----------
    coords :
        Coordinates of the source
    f :
        Source function (Ricker)

    domainFlag :
        0 the source is not in the domain,
        1 the source is in the domain
    """
    def __init__(self,
                 coords=None,
                 source=None):
        """Constructor for the point source

        Parameters
        ----------
        coords :
            Coordinates for the point source

        source :
            Source function (Ricker)
        """

        self.coords = coords
        self.value = source

    def __repr__(self):
        return '('+str(self.coords[0])+','+str(self.coords[1])+','+str(self.coords[2])+')'


    def x(self):
        return self.coords[0]

    def y(self):
        return self.coords[1]

    def z(self):
        return self.coords[2]


     def construct_from_dict(self, **kwargs):
        for key, value in kwargs.items():
            setattr(self, key, value)

        return self


class SourceSet:
    """A class representing a point source

    Attributes
    ----------
    source_list :
        List of Source object
    n :
        Number of Source
    """

    def __init__(self,
                 source_list=None):

        if source_list is None:
            self.source_list = []
            self.n = 0
        else:
            self.source_list = source_list
            self.n = len(source_list) #Number of sources


    def append(self, source):
        self.source_list.append(source)

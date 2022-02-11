from copy import deepcopy
import numpy as np
#import segyio
import random
import os
import xml.etree.ElementTree as ET

class Acquisition:
    def __init__(self,
                 dt=None,
                 velocity_model=None,
                 acquisition_type=None,
                 boundary=None,
                 source=None,
                 source_positions=None,
                 receivers_positions=None,
                 source_depth=None,
                 receivers_depth=None,
                 output=None,
                 **kwargs):

        self.xml=None
        self.limited_aperture=False
        self.boundary=boundary
        self.velocity_model=velocity_model
        self.aperture_boundary=boundary

        if source_positions is not None and receivers_positions is not None:
            self.acquisition(source,
                             source_positions,
                             receivers_positions,
                             source_depth,
                             receivers_depth)

            for shot in self.shots:
                shot.velocity_model = velocity_model

        if dt is not None:
            for shot in self.shots:
                shot.dt = dt

        if output is not None:
            self.output = output
        else:
            self.output = os.path.join(os.getcwd(),"seismoTrace/")
            if os.path.exists(self.output):
                pass
            else:
                os.system("mkdir -p " + self.output)


    def __repr__(self):
        string_list = []
        string_list.append("Acquisition type : " + type(self).__name__ + "\n")
        string_list.append("Number of shots : " + str(len(self.shots)) + "\n")
        string_list.append("Acquisition boundaries : " + str(self.boundary) +"\n")
        if self.limited_aperture == True:
            string_list.append("Limited aperture \n")
        else:
            string_list.append("Full aperture \n")
        if isinstance(self.velocity_model, int):
            string_list.append("Homogeneous velocity model \n")
        else:
            string_list.append("Heterogeneous velocity model \n")

        rep=""
        for string in string_list:
            rep += string

        return rep


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
        nb_shot_m1 = len(self.shots)
        ind = 0

        for i in range(n):
            nb_shot = int(nb_shot_m1/(n-i))

            acqs.append(deepcopy(self))
            acqs[i].shots = acqs[i].shots[ind:ind + nb_shot]

            ind = ind + nb_shot
            nb_shot_m1 = nb_shot_m1 - nb_shot

        return acqs



    def setDt(self, dt):
        for shot in self.shots:
            shot.dt=dt


    #Can't use it for now
    def calculDt(self):
        if self.xml is None:
            raise ValueError("You must link the seismic acquisition with xml file to be able to calculte the simulation time step dt")
        else:
            for shot in self.shots:
                tree = ET.parse(shot.xml)
                root = tree.getroot()

                internal_mesh = [elem.attrib for elem in root.iter('InternalMesh')][0]

                xCoords_str_list = internal_mesh['xCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()
                yCoords_str_list = internal_mesh['yCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()
                zCoords_str_list = internal_mesh['zCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()

                xCoords = [float(xCoords_str_list[0]), float(xCoords_str_list[1])]
                yCoords = [float(yCoords_str_list[0]), float(yCoords_str_list[1])]
                zCoords = [float(zCoords_str_list[0]), float(zCoords_str_list[1])]

                x_nb_cells = int(internal_mesh['nx'].replace('{','').replace('}',''))
                y_nb_cells = int(internal_mesh['ny'].replace('{','').replace('}',''))
                z_nb_cells = int(internal_mesh['nz'].replace('{','').replace('}',''))

                x = (xCoords[1] - xCoords[0])/x_nb_cells
                y = (yCoords[1] - yCoords[0])/y_nb_cells
                z = (zCoords[1] - zCoords[0])/z_nb_cells

                minRadius = min(x/2, y/2, z/2)

                if isinstance(shot.velocity_model, str):

                    with open(shot.velocity_model,'r') as vf:
                        minVelocity = float(vf.readline())
                        for velocity in vf.readlines():
                            if float(velocity) < minVelocity:
                                minVelocity = float(velocity)
                    vf.close()

                else:
                    minVelocity = shot.velocity_model

                #Suppose space discretization is order 1
                shot.dt = minRadius/minVelocity



    def limitedAperture(self, aperture_dist):
        self.limited_aperture=True

        tree = ET.parse(self.xml)
        root = tree.getroot()

        internal_mesh = [elem.attrib for elem in root.iter('InternalMesh')][0]
        box = [elem.attrib for elem in root.iter('Box')][0]
        if len([elem.attrib for elem in root.iter('VTK')]):
            vtk = [elem.attrib for elem in root.iter('VTK')][0]
        else:
            vtk = []
        for periodicEvent in [elem.attrib for elem in root.iter('PeriodicEvent')]:
            if periodicEvent['name'] == "vtk":
                events = periodicEvent
            else:
                events = []

        xCoords_str_list = internal_mesh['xCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()
        yCoords_str_list = internal_mesh['yCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()
        zCoords_str_list = internal_mesh['zCoords'].strip('').replace(',', ' ').replace('{','').replace('}','').split()

        xCoords = [float(xCoords_str_list[0]), float(xCoords_str_list[1])]
        yCoords = [float(yCoords_str_list[0]), float(yCoords_str_list[1])]
        zCoords = [float(zCoords_str_list[0]), float(zCoords_str_list[1])]

        x_nb_cells = int(internal_mesh['nx'].replace('{','').replace('}',''))
        y_nb_cells = int(internal_mesh['ny'].replace('{','').replace('}',''))
        z_nb_cells = int(internal_mesh['nz'].replace('{','').replace('}',''))

        x_cells_boundary = np.append(np.arange(xCoords[0], xCoords[1], (xCoords[1]-xCoords[0])/x_nb_cells), xCoords[1])
        y_cells_boundary = np.append(np.arange(yCoords[0], yCoords[1], (yCoords[1]-yCoords[0])/y_nb_cells), yCoords[1])
        z_cells_boundary = np.append(np.arange(zCoords[0], zCoords[1], (zCoords[1]-zCoords[0])/z_nb_cells), zCoords[1])

        self.set_limited_aperture_boundaries(aperture_dist,
                                             x_cells_boundary,
                                             y_cells_boundary,
                                             z_cells_boundary)

        if os.path.exists("limitedAperture") == False:
            os.system("mkdir -p limitedAperture")


        if isinstance(self.velocity_model, str):
            tablefunction =[elem.attrib for elem in root.iter('TableFunction')][0]

            with open(self.velocity_model.split("_velModel")[0]+"_xlin.geos", 'r') as xf:
                xlin = xf.readlines()
            xf.close()
            lx = len(xlin)
            with open(self.velocity_model.split("_velModel")[0]+"_ylin.geos", 'r') as yf:
                ylin = yf.readlines()
            yf.close()
            ly = len(ylin)
            with open(self.velocity_model.split("_velModel")[0]+"_zlin.geos", 'r') as zf:
                zlin = zf.readlines()
            zf.close()
            lz = len(zlin)

            velocity_model_value = []
            velocity_model_index = []
            with open(self.velocity_model,'r') as vf:
                for ind, line in enumerate(vf):
                    velocity_model_value.append(line)
                    velocity_model_index.append(ind)
            vf.close()


        for shot in self.shots:
            nx = int(x_nb_cells*(shot.boundary[0][1]-shot.boundary[0][0])/(xCoords[1]-xCoords[0]))
            ny = int(y_nb_cells*(shot.boundary[1][1]-shot.boundary[1][0])/(yCoords[1]-yCoords[0]))
            nz = int(z_nb_cells*(shot.boundary[2][1]-shot.boundary[2][0])/(zCoords[1]-zCoords[0]))

            internal_mesh['xCoords'] = "{"+str(shot.boundary[0][0])+","+str(shot.boundary[0][1])+"}"
            internal_mesh['yCoords'] = "{"+str(shot.boundary[1][0])+","+str(shot.boundary[1][1])+"}"
            internal_mesh['zCoords'] = "{"+str(shot.boundary[2][0])+","+str(shot.boundary[2][1])+"}"

            internal_mesh['nx'] = "{"+str(nx)+"}"
            internal_mesh['ny'] = "{"+str(ny)+"}"
            internal_mesh['nz'] = "{"+str(nz)+"}"

            box['xMin'] = "{"+str(shot.boundary[0][0]-0.01)+","+str(shot.boundary[1][0]-0.01)+","+str(shot.boundary[2][1]-0.01)+"}"
            box['xMax'] = "{"+str(shot.boundary[0][1]+0.01)+","+str(shot.boundary[1][1]+0.01)+","+str(shot.boundary[2][1]+0.01)+"}"

            if isinstance(self.velocity_model, str):
                xlocal=[]
                ylocal=[]
                zlocal=[]

                keep_line=[]
                localToGlobal=[]
                for i in range(lx):
                    if float(xlin[i])>shot.boundary[0][0] or float(xlin[i])<shot.boundary[0][1]:
                        xlocal.append(xlin[i])
                        for j in range(ly):
                            if float(ylin[j])>shot.boundary[1][0] or float(ylin[j])<shot.boundary[1][1]:
                                if ylin[j] not in ylocal:
                                    ylocal.append(ylin[j])
                                for k in range(lz):
                                    if float(zlin[k])>shot.boundary[2][0] or float(zlin[k])<shot.boundary[2][1]:
                                        if zlin[k] not in zlocal:
                                            zlocal.append(zlin[k])
                                        keep_line.append(velocity_model_value[i*ly*lz+j*lz+k])
                                        localToGlobal.append(i*ly*lz+j*lz+k)

                xlin_file = self.velocity_model.split("_velModel")[0]+"_xlin"+shot.id+".geos"
                ylin_file = self.velocity_model.split("_velModel")[0]+"_ylin"+shot.id+".geos"
                zlin_file = self.velocity_model.split("_velModel")[0]+"_zlin"+shot.id+".geos"

                with open(xlin_file,'w') as xfla:
                    for x in xlocal:
                        xfla.write(x)
                xfla.close()
                with open(ylin_file,'w') as yfla:
                    for y in ylocal:
                        yfla.write(y)
                yfla.close()
                with open(zlin_file,'w') as zfla:
                    for z in zlocal:
                        zfla.write(z)
                zfla.close()

                velocity_file = self.velocity_model.split('.')[0] + shot.id+ ".geos"

                max_vel=0
                with open(velocity_file, 'w') as vfla:
                    for line in keep_line:
                        vfla.write(line)
                        if float(line) > max_vel:
                            max_vel =  float(line)
                vfla.close()

                tablefunction['coordinateFiles']="{"+xlin_file+","+ylin_file+","+zlin_file+"}"
                tablefunction['voxelFile'] = velocity_file
                shot.velocity_model = velocity_file

            if len(vtk) != 0:
                vtk['name'] = 'vtkOutput'+str(shot.id)
            if len(events) != 0:
                events['target'] = '/Outputs/vtkOutput'+str(shot.id)

            xmlfile = os.path.join("limitedAperture/", "limited_aperture_Shot"+str(shot.id)+".xml")
            tree.write(xmlfile)

            shot.xml = xmlfile



    def set_limited_aperture_boundaries(self,
                                        distance,
                                        x_cells_boundary,
                                        y_cells_boundary,
                                        z_cells_boundary):

        for shot in self.shots:
            nx = x_cells_boundary.size
            xmin = x_cells_boundary[0]
            xmax = x_cells_boundary[nx-1]

            for i in range(nx):
                if shot.sources.source_list[0].x() - distance >= x_cells_boundary[i]:
                    xmin = x_cells_boundary[i]
                if shot.sources.source_list[0].x() + distance <= x_cells_boundary[nx-1-i]:
                    xmax = x_cells_boundary[nx-1-i]

            ny = y_cells_boundary.size
            ymin = y_cells_boundary[0]
            ymax = y_cells_boundary[ny-1]

            for i in range(ny):
                if shot.sources.source_list[0].y() - distance >= y_cells_boundary[i]:
                    ymin = y_cells_boundary[i]
                if shot.sources.source_list[0].y() + distance <= y_cells_boundary[ny-1-i]:
                    ymax = y_cells_boundary[ny-1-i]

            nz = z_cells_boundary.size
            zmin = z_cells_boundary[0]
            zmax = z_cells_boundary[nz-1]

            limited_aperture = [[xmin,xmax],[ymin,ymax],[zmin,zmax]]

            shot.boundary = deepcopy(limited_aperture)
            shot.receivers.keep_inside_domain(limited_aperture)


    def add_xml(self, xmlfile):
        self.xml = xmlfile
        for shot in self.shots:
            shot.xml = xmlfile



    def construct_from_dict(self, **kwargs):
        self.shots=[]
        for key, value in kwargs.items():
            if key == "shots":
                for i in range(len(value)):
                    dict = value[i]
                    self.shots.append(Shot().construct_from_dict(**dict))
            else:
                setattr(self, key, value)

        return self



"""
class SEGYAcquisition(Acquisition):

    def __init__(self,
                 limited_aperture=False,
                 velocity_model=None,
                 boundary=None,
                 dt=None,
                 directory=None,
                 source=None):

        if directory is not None:
            self.acquisition(directory,
                             source)

            for shot in self.shots:
                shot.velocity_model = velocity_model
        else:
            raise ValueError("You must specify a directory")

        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model,
                         dt=dt)



    def acquisition(self,
                    directory,
                    wavelet):

        Acquisition from .sgy files

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

                shots.append(Shot(source, receivers, shot_id))
                j+=1
            i+=1

        self.shots = shots

"""

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

            for shot in self.shots:
                shot.velocity_model = velocity_model


        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model,
                         dt=dt)



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
        receivers.append(deepcopy(receivers_temp))

        xs = np.linspace(start_source_pos[0], end_source_pos[0], number_of_sources).tolist()
        ys = np.linspace(start_source_pos[1], end_source_pos[1], number_of_sources).tolist()

        shots = []

        for i in range(len(xs)):
            source = SourceSet()
            if len(str(i+1))<2:
                shot_id = "00"+str(i+1)
            elif len(str(i+1))<3 and len(str(i+1))>2:
                shot_id = "0"+str(i+1)
            else:
                shot_id = str(i+1)

            srcpos = [xs[i], ys[i], source_depth]
            source.append(Source(srcpos, wavelet))

            shot = Shot(source, receivers, shot_id)
            shots.append(deepcopy(shot))

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

            for shot in self.shots:
                shot.velocity_model = velocity_model


        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         source=source,
                         velocity_model=velocity_model,
                         dt=dt)



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
                shots.append(deepcopy(shot))


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


        if number_of_sources is not None and number_of_receivers_x is not None and number_of_receivers_y is not None:
            self.acquisition(wavelet=source,
                             boundary=boundary,
                             number_of_sources=number_of_sources,
                             number_of_receivers_x=number_of_receivers_x,
                             number_of_receivers_y=number_of_receivers_y,
                             source_depth=source_depth,
                             receivers_depth=receivers_depth)

            for shot in self.shots:
                shot.velocity_model = velocity_model


        super().__init__(limited_aperture=limited_aperture,
                         boundary=boundary,
                         velocity_model=velocity_model,
                         dt=dt,
                         source=source)



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
            shots.append(deepcopy(shot))

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
                 sources=None,
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

        self.sources = sources
        self.receivers = receivers
        self.flag = "Undone"
        self.dt = None
        self.xml = None
        self.velocity_model = None
        self.boundary = None

        if shot_id is not None:
            self.id = shot_id
        else:
            self.id = None

    def __repr__(self):
        return 'Source position : \n'+str(self.sources) +' \n\n' + 'Receivers positions : \n' + str(self.receivers) + '\n\n'


    def flagUpdate(self, string):
        self.flag = string


    def construct_from_dict(self, **kwargs):
        for key, value in kwargs.items():
            if key == "source":
                setattr(self, key, Source().construct_from_dict(**value))
            elif key == "receivers":
                setattr(self, key, ReceiverSet().construct_from_dict(**value))
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

    def __init__(self, receivers_list = None):
        """Constructor for the point source

        Parameters
        ----------
        receiver_list :
            List of Receiver
        """
        if receivers_list is None:
            self.receivers_list = []
            self.n = 0
        else:
            self.receivers_list = receivers_list
            self.n = len(receivers_list)

    def __repr__(self):
        if self.n >=10:
            return str(self.receivers_list[0:4])[:-1] + '...' + '\n' + str(self.receivers_list[-4:])[1:]
        else:
            return str(self.receivers_list)


    def keep_inside_domain(self, boundaries):
        new_list=[]

        for receiver in self.receivers_list:
            if boundaries[0][0] <= receiver.x() and receiver.x() <= boundaries[0][1] \
            and boundaries[1][0] <= receiver.y() and receiver.y() <= boundaries[1][1] \
            and boundaries[2][0] <= receiver.z() and receiver.z() <= boundaries[2][1]:
                new_list.append(receiver)

        self.receivers_list = new_list
        self.n = len(new_list)


    def getSetCoord(self):
        rcvs_coords_list = []
        for i in range(self.n):
            rcvs_coords_list.append(np.array([self.receivers_list[i].coords[0], self.receivers_list[i].coords[1], self.receivers_list[i].coords[2]]))
        return rcvs_coords_list

    def append(self, ReceiverSet):
        for i in range(ReceiverSet.n):
            self.receivers_list.append(ReceiverSet.receivers_list[i])
        self.n += ReceiverSet.n


    #Transpose all receivers position
    def linearSetMove(self, coord, value):
        for i in range(self.n):
            self.receivers_list[i].linearMove(coord, value)


    def construct_from_dict(self, **kwargs):
        self.receivers_list = []
        for key, value in kwargs.items():
            if key == "receivers_list":
                for i in range(len(value)):
                    dict = value[i]
                    self.receivers_list.append(Receiver().construct_from_dict(**dict))
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
            self.n = len(source_list)

    def __repr__(self):
        if self.n >=10:
            return str(self.source_list[0:4])[:-1] + '...' + '\n' + str(self.source_list[-4:])[1:]
        else:
            return str(self.source_list)


    def append(self, source):
        self.source_list.append(source)

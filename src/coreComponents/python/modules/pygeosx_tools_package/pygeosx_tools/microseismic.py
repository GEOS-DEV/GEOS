import numpy as np
from pygeosx_tools import wrapper, parallel_io
from scipy.spatial import Delaunay
from pyevtk.hl import pointsToVTK


class JointSet():
    """
    Joint used for microseismic analysis.
    Note: values that are in the strike-dip-normal reference
          frame are indicated by 'sdn'

    Attributes:
      region_name (str): GEOSX mesh region name
      solid_name (str): GEOSX solid material name
      fluid_name (str): GEOSX fluid material name
      sigma_key (str): Key pointing to the GEOSX stress wrapper
      pressure_key (str): Key pointing to the GEOSX pressure wrapper
      allow_repeating (bool): Flag to allow repeating events on joints
      shear_modulus (float): Shear modulus of host rock (Pa)
      N (int): Number of joints in the set
      strike (list): Strike (counter-clockwise from x-axis, radians)
      dip (list): Dip (radians)
      mu (list): Coefficient of friction (unitless)
      cohesion (list): Cohesion (Pa)
      element_id (list): Parent element index
      location (list): Joint locations (m)
      sigma_sdn_offset (list): Stress offset in sdn frame (Pa)
      fail_count (list): Cumulative number of microseismic events on the joint
      rotation (list): Rotation matrices to map from xyz to sdn frames
    """

    def __init__(self,
                 region_name='Region',
                 solid_name='rock',
                 fluid_name='',
                 allow_repeating=False,
                 shear_modulus=1.0e9):
        # GEOSX object names / keys
        self.region_name = region_name
        self.solid_name = solid_name
        self.fluid_name = fluid_name
        self.sigma_key = ''
        self.pressure_key = ''

        # Joint behavior
        self.allow_repeating = allow_repeating
        self.shear_modulus = shear_modulus

        # Joint values
        self.N = 0
        self.strike = []
        self.dip = []
        self.mu = []
        self.cohesion = []
        self.element_id = []
        self.location = []
        self.sigma_sdn_offset = []
        self.fail_count = []
        self.rotation = []

    def build_joint_set(self,
                        problem,
                        random_location=True,
                        density=1.0,
                        uniform_orientation=True,
                        strike_mean=0.0,
                        strike_std=0.0,
                        dip_mean=0.0,
                        dip_std=0.0,
                        mu_mean=0.6,
                        mu_std=0.0,
                        cohesion_mean=1.0e6,
                        cohesion_std=0.0):
        """
        Add a joint set to a given region

        Args:
            problem (pygeosx.group): GEOSX problem handle
            random_location (bool): Choose random locations within the regin (default=True)
            density (float): Joint density (#/m3)
            uniform_orientation (bool): Use randomly oriented joints (default=True)
            strike_mean (float): Strike angle mean (degrees)
            strike_std (float): Strike angle std deviation (degrees)
            dip_mean (float): Dip angle mean (degrees)
            dip_std: Dip angle std deviation (degrees)
            mu_mean (float): Friction coefficient mean (unitless)
            mu_std: Friction coefficient std deviation (unitless)
            cohesion_mean (float): Cohesion mean (Pa)
            cohesion_std: Cohesion std deviation (Pa)
        """
        # Choose locations
        if random_location:
            self.choose_random_locations(problem, density)
        else:
            # TODO: Add structured joint set
            raise Exception('Option not available: random_location=False')

        # Choose orientations
        if uniform_orientation:
            self.choose_random_uniform_orientations()
        else:
            self.choose_random_normal_orientations(strike_mean=strike_mean,
                                                   strike_std=strike_std,
                                                   dip_mean=dip_mean,
                                                   dip_std=dip_std)

        # Setup other values
        self.mu = (np.random.randn(self.N) * mu_std) + mu_mean
        self.cohesion = (np.random.randn(self.N) * cohesion_std) + cohesion_mean
        self.rotation = [self.get_rotation_matrix(s, d) for s, d in zip(self.strike, self.dip)]
        self.sigma_sdn_offset = np.zeros((self.N, 3, 3))
        self.fail_count = np.zeros(self.N, dtype=int)
        self.sigma_key = wrapper.get_matching_wrapper_path(problem, [self.region_name, self.solid_name, 'stress'])

        if self.fluid_name:
            self.pressure_key = wrapper.get_matching_wrapper_path(problem, [self.region_name, self.fluid_name, 'pressure'])

    def choose_random_locations(self, problem, density):
        """
        Choose random joint locations within the target region

        Args:
            problem (pygeosx.group): GEOSX problem handle
            density (float): Joint density (#/m3)
        """
        # Find key values
        ghost_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'ghostRank'])
        node_position_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'ReferencePosition'])
        node_element_key = wrapper.get_matching_wrapper_path(problem, ['nodeManager', 'elemList'])
        element_node_key = wrapper.get_matching_wrapper_path(problem, [self.region_name, 'nodeList'])

        # Find the region extents
        g = wrapper.get_wrapper(problem, ghost_key)
        x = wrapper.get_wrapper(problem, node_position_key)
        xb = x[g < 0, :]
        xmin = np.amin(xb, axis=0)
        xmax = np.amax(xb, axis=0)

        # Choose the inital number of elements to test connectivity
        Ntest = int(np.prod(xmax - xmin) * density)
        tmp_index = np.zeros(Ntest, dtype=int)
        tmp_location = np.zeros((Ntest, 3))
        test_locations = np.random.uniform(size=(Ntest, 3))
        for ii in range(3):
            test_locations[:, ii] *= (xmax[ii] - xmin[ii])
            test_locations[:, ii] += xmin[ii]

        # Map joint location to elements
        node_element_map = wrapper.get_wrapper(problem, node_element_key)
        element_node_map = wrapper.get_wrapper(problem, element_node_key)
        self.N = 0
        for ii in range(Ntest):
            # Find the closest node
            R2 = (x[:, 0] - test_locations[ii, 0])**2 + (x[:, 1] - test_locations[ii, 1])**2 + (x[:, 2] - test_locations[ii, 2])**2
            Ia = np.argmin(R2)

            # Test to see if the joint is within an adjoining element
            test_elements = node_element_map[Ia]
            for element in test_elements:
                test_nodes = element_node_map[element]
                element_corners = x[test_nodes, :]
                if self.point_in_convex_element(test_locations[ii, :], element_corners):
                    tmp_index[self.N] = element
                    tmp_location[self.N, :] = test_locations[ii, :]
                    self.N += 1
                    break

        # Finalize location selection
        print('%i joints generated' % (self.N))
        self.element_id = tmp_index[:self.N]
        self.location = tmp_location[:self.N, :]

    def point_in_convex_element(self, point, vertices):
        """
        Check whether a point is within a convex element

        Args:
            point (np.array): target location (N)
            vertices (np.array): element vertices (MxN)
        """
        hull = Delaunay(vertices)
        return hull.find_simplex(point) >= 0

    def check_failure_criteria(self, problem):
        """
        Check the failure criteria on each joint

        Args:
            problem (pygeosx.group): GEOSX problem handle
        """
        events = []
        sigma = wrapper.get_wrapper(problem, self.sigma_key)
        pressure = []
        if self.pressure_key:
            pressure = wrapper.get_wrapper(problem, self.pressure_key)
        for ii in range(self.N):
            s = np.zeros((3, 3))
            if ((self.fail_count[ii] == 0) | self.allow_repeating):
                # Rotate to the sdn frame
                jj = self.element_id[ii]
                s[0, 0] = sigma[jj, 0, 0]
                s[1, 1] = sigma[jj, 0, 1]
                s[2, 2] = sigma[jj, 0, 2]
                s[0, 1] = sigma[jj, 0, 3]
                s[0, 2] = sigma[jj, 0, 4]
                s[1, 2] = sigma[jj, 0, 5]
                s[1, 0] = s[0, 1]
                s[2, 0] = s[0, 2]
                s[2, 1] = s[1, 2]
                R = self.rotation[ii]
                s_sdn = np.dot(np.dot(R, s), np.transpose(R))
                s_sdn -= self.sigma_sdn_offset[ii, ...]

                # Check failure criteria
                tau = np.sqrt(s_sdn[0, 2]**2 + s_sdn[1, 2]**2)
                sigma_effective = s_sdn[2, 2]
                if pressure:
                    sigma_effective -= pressure[ii]

                if (tau > (self.cohesion[ii] - self.mu[ii] * sigma_effective)):
                    self.fail_count[ii] += 1
                    rake = np.arctan2(s_sdn[1, 2], s_sdn[0, 2])
                    dT = 0.001 * self.shear_modulus
                    self.sigma_sdn_offset[ii, 0, 2] += dT * np.cos(rake)
                    self.sigma_sdn_offset[ii, 1, 2] += dT * np.sin(rake)
                    events.append([self.location[ii, ...],
                                   self.strike[ii],
                                   self.dip[ii],
                                   rake])
        return events

    def choose_random_uniform_orientations(self):
        """
        Choose random orientations for joints that
        are uniform across the focal sphere
        """
        self.strike = np.random.uniform(low=0.0,
                                        high=2.0*np.pi,
                                        size=self.N)
        self.dip = np.arccos(np.random.uniform(low=-1.0,
                                               high=1.0,
                                               size=self.N))

    def choose_random_normal_orientations(self, strike_mean, strike_std, dip_mean, dip_std):
        """
        Choose random normal orientations for joints

        Args:
            strike_mean (float): Strike angle mean (degrees)
            strike_std (float): Strike angle std deviation (degrees)
            dip_mean (float): Dip angle mean (degrees)
            dip_std: Dip angle std deviation (degrees)
        """
        self.strike = (np.random.randn(self.N) * strike_std) + strike_mean
        self.dip = (np.random.randn(self.N) * dip_std) + dip_mean

    def get_rotation_matrix(self, strike, dip):
        """
        Get the rotation matrix to convert from xyz to sdn frames

        Args:
            strike (float): Strike (counter-clockwise from x-axis, radians)
            dip (float): Dip (radians)

        Returns:
            np.array: Rotation matrix
        """
        stheta = np.sin(strike * np.pi / 180.0)
        ctheta = np.cos(strike * np.pi / 180.0)
        sdelta = np.sin(dip * np.pi / 180.0)
        cdelta = np.cos(dip * np.pi / 180.0)

        rotation_matrix = np.zeros((3, 3))
        rotation_matrix[0, 0] = ctheta
        rotation_matrix[1, 0] = stheta
        rotation_matrix[2, 0] = 0.0
        rotation_matrix[0, 1] = -stheta*cdelta
        rotation_matrix[1, 1] = ctheta*cdelta
        rotation_matrix[2, 1] = -sdelta
        rotation_matrix[0, 2] = -stheta*sdelta
        rotation_matrix[1, 2] = ctheta*sdelta
        rotation_matrix[2, 2] = cdelta

        return rotation_matrix


class MicroseismicAnalysis():
    def __init__(self):
        """
        Microseismic Analysis class
        """
        self.joint_sets = []
        self.set_names = []
        self.catalog_order = ['t',
                              'x',
                              'y',
                              'z',
                              'strike',
                              'dip',
                              'rake',
                              'set_id']
        self.catalog = {k: [] for k in self.catalog_order}
        self.catalog_requires_sync = False
        self.catalog_all = {}

    def add_joint_set(self, problem, set_name, **xargs):
        """
        Add a joint set to a given region

        Args:
            problem (pygeosx.group): GEOSX problem handle
            set_name (str): Joint set name
            xargs: See JointSet.__init__ and JointSet.build_joint_set
        """
        init_args = ['region_name', 'solid_name', 'fluid_name', 'allow_repeating', 'shear_modulus']
        joint_xargs = {k: xargs[k] for k in xargs.keys() if k in init_args}
        build_xargs = {k: xargs[k] for k in xargs.keys() if k not in init_args}

        self.joint_sets.append(JointSet(**joint_xargs))
        self.set_names.append(set_name)
        self.joint_sets[-1].build_joint_set(problem, **build_xargs)

    def check_microseismic_activity(self, problem):
        """
        Check microseismic activity for each joint set,
        and save any events to the catalog

        Args:
            problem (pygeosx.group): GEOSX problem handle
        """
        time = wrapper.get_wrapper(problem, 'Events/time')[0]
        self.catalog_requires_sync = True
        for ii, joint_set in enumerate(self.joint_sets):
            events = joint_set.check_failure_criteria(problem)
            for event in events:
                self.catalog['t'].append(time)
                self.catalog['x'].append(event[0][0])
                self.catalog['y'].append(event[0][1])
                self.catalog['z'].append(event[0][2])
                self.catalog['strike'].append(event[1] * 180.0 / np.pi)
                self.catalog['dip'].append(event[2] * 180.0 / np.pi)
                self.catalog['rake'].append(event[3] * 180.0 / np.pi)
                self.catalog['set_id'].append(ii)

    def sync_catalog(self):
        if self.catalog_requires_sync:
            self.catalog_requires_sync = False
            if (parallel_io.comm.size == 1):
                # This is a serial run
                self.catalog_all = {k: np.ascontiguousarray(v) for k, v in self.catalog.items()}
            else:
                for k, v in self.catalog.items():
                    self.catalog_all[k] = np.ascontiguousarray(parallel_io.gather_array(v))

    def save_catalog_csv(self, fname):
        """
        Save the catalog to a csv format file

        Args:
            fname (str): Name of the file (with no extension)
        """
        self.sync_catalog()
        if (parallel_io.rank == 0):
            catalog_header = ', '.join(self.catalog_order)
            catalog_data = np.transpose(np.array([self.catalog_all[k] for k in self.catalog_order]))
            np.savetxt('%s.csv' % (fname),
                       catalog_data,
                       fmt='%1.4f',
                       header=catalog_header,
                       delimiter=', ')

    def save_catalog_vtk(self, fname):
        """
        Save the catalog to a vtk format file

        Args:
            fname (str): Name of the file (with no extension)
        """
        self.sync_catalog()
        if (parallel_io.rank == 0):
            pointsToVTK(fname,
                        self.catalog_all['x'],
                        self.catalog_all['y'],
                        self.catalog_all['z'],
                        data=self.catalog_all)



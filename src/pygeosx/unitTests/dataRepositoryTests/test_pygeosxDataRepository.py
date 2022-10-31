import pylvarray
import pygeosx
import numpy as np
import pytest
import os

# Geosx paths
level_zero = 'domain/MeshBodies/mesh1/meshLevels/Level0/'
node_manager = level_zero + 'nodeManager/'
cell_blocks = level_zero + 'ElementRegions/elementRegionsGroup/dummy/elementSubRegions/cellBlock01/'


class TestPygeosxObjects(object):

    def get_geosx_args(self):
        script_path = os.path.dirname(os.path.abspath(__file__))
        target_xml = os.path.join(script_path, 'pygeosx_only.xml')
        return ('pygeosx', '-i', target_xml)

    @pytest.fixture(scope='class')
    def geosx_problem(self):
        """
        Initialize the minimal pygeosx problem that has the following features:
            A single rank
            A single periodic event with dt=1 targeting a python output event

        Returns:
            pygeosx.problem: the problem handle
        """
        args = self.get_geosx_args()
        problem = pygeosx.initialize(0, args)
        return problem

    def test_apply_initial_conditions(self, geosx_problem):
        pygeosx.apply_initial_conditions()

    def test_get_group(self, geosx_problem):
        """
        Test accessing a group, which has an expected wrapper 'time'
        """
        group = geosx_problem.get_group('Events', None)
        assert isinstance(group, pygeosx.Group)
        t = group.get_wrapper('time').value()[0]
        assert isinstance(t, float)

    # yapf: disable
    @pytest.mark.parametrize('target_key, target_shape, target_type',
                             [('Events/time', (1, ), float),
                              ('Events/cycle', (1, ), np.int32),
                              ('Events/python/target', 21, str),
                              (node_manager + 'ReferencePosition', (8, 3), float),
                              (node_manager + 'ghostRank', (8, ), np.int32),
                              (cell_blocks + 'elementVolume', (1, ), float)])
    def test_wrapper_read(self, target_key, target_shape, target_type, geosx_problem):
        """
        Test whether wrappers can be read, and whether their shape and value types are correct

        Args:
            target_key (str): Key value for the target wrapper
            target_shape (tuple): Expected shape for the target wrapper
            target_type (dtype): Expected type for the wrapper elements
            geosx_proplem (pygeosx.problem): The problem handle
        """
        # yapf: enable
        s = 0
        val = 0
        x = geosx_problem.get_wrapper(target_key).value()
        if hasattr(x, "to_numpy"):
            x = x.to_numpy()

        if isinstance(x, np.ndarray):
            s = np.shape(x)
            Ia = tuple([0 for y in s])
            val = x[Ia]
        else:
            s = len(x)
            val = x

        assert s == target_shape
        assert isinstance(val, target_type)

    # yapf: disable
    @pytest.mark.parametrize('target_key',
                             [('Events/time'),
                              ('Events/cycle'),
                              (node_manager + 'ReferencePosition'),
                              (node_manager + 'ghostRank'),
                              (cell_blocks + 'elementVolume')])
    def test_wrapper_write(self, target_key, geosx_problem):
        """
        Test writing to the wrappers

        Args:
            target_key (str): Key value for the target wrapper
            geosx_proplem (pygeosx.problem): The problem handle
        """
        # yapf: enable
        # Access and record the initial value
        x = geosx_problem.get_wrapper(target_key).value()
        if hasattr(x, "set_access_level"):
            x.set_access_level(pylvarray.MODIFIABLE, pylvarray.CPU)
        if hasattr(x, "to_numpy"):
            x = x.to_numpy()
        x_init = x.copy()

        # Modify the value
        x += 1
        del x

        # Get another copy of the wrapper
        x = geosx_problem.get_wrapper(target_key).value()
        if hasattr(x, "to_numpy"):
            x = x.to_numpy()

        # Check to see if the value changed appropriately
        dx = x_init + 1 - x
        if isinstance(dx, np.ndarray):
            dx = np.sum(dx)
        dx = float(dx)

        assert dx < 1e-10

    # yapf: disable
    @pytest.mark.parametrize('target_key',
                             [(node_manager + 'ReferencePosition'),
                              (node_manager + 'ghostRank'),
                              (cell_blocks + 'elementVolume')])
    def test_wrapper_read_only_error(self, target_key, geosx_problem):
        """
        Test writing to read-only wrappers.
        Note: this error does not get called on certain scalar wrappers (time, cycle)

        Args:
            target_key (str): Key value for the target wrapper
            geosx_proplem (pygeosx.problem): The problem handle
        """
        # yapf: enable
        x = geosx_problem.get_wrapper(target_key).value()
        if hasattr(x, "to_numpy"):
            print('Test_a')
            x = x.to_numpy()

        with pytest.raises(ValueError):
            x += 1

    def test_problem_run(self, geosx_problem):
        """
        Test advancing the problem one step
        """
        status = pygeosx.run()
        assert status == pygeosx.READY_TO_RUN
        group = geosx_problem.get_group('Events', None)
        ta = group.get_wrapper('time').value()[0].copy()
        ca = group.get_wrapper('cycle').value()[0].copy()

        status = pygeosx.run()
        assert status == pygeosx.READY_TO_RUN
        tb = group.get_wrapper('time').value()[0]
        cb = group.get_wrapper('cycle').value()[0]
        assert (ta + 1) == tb
        assert (ca + 1) == cb

    def test_problem_complete(self, geosx_problem):
        """
        Test running the problem to completion
        """
        max_iter = 100
        ii = 0
        while (pygeosx.run() != pygeosx.COMPLETED) and (ii <= max_iter):
            ii += 1

        assert ii != max_iter

    def test_problem_reinit(self, geosx_problem):
        """
        Test re-initializing the problem
        """
        args = self.get_geosx_args()
        problem = pygeosx.reinit(args)
        c = problem.get_wrapper('Events/cycle').value()[0]
        assert c == 0

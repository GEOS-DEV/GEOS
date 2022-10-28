
import pylvarray
import pygeosx
import numpy as np
import pytest
import os


class TestPygeosxObjects(object):

    @pytest.fixture(scope='class')
    def geosx_problem(self):
        """
        Initialize the minimal pygeosx problem that has the following features:
            A single rank
            A single periodic event with dt=1 targeting a python output event

        Returns:
            pygeosx.problem: the problem handle
        """
        script_path = os.path.dirname(os.path.abspath(__file__))
        target_xml = os.path.join(script_path, 'pygeosx_only.xml')

        rank = 0
        args = ('pygeosx', '-i', target_xml)
        problem = pygeosx.initialize(rank, args)
        pygeosx.apply_initial_conditions()
        # print(wrapper.get_matching_wrapper_path(problem, ['cellBlock01']))
        return problem

    def test_get_group(self, geosx_problem):
        """
        Test accessing a group, which has an expected wrapper 'time'
        """
        group = geosx_problem.get_group('Events', None)
        t = group.get_wrapper('time')

    @pytest.mark.parametrize('target_key, target_shape, target_type',
                             [('Events/time', (1,), float),
                              ('Events/cycle', (1,), np.int32),
                              ('Events/python/target', 21, str),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ReferencePosition', (8, 3), float),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ghostRank', (8,), np.int32),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/ElementRegions/elementRegionsGroup/dummy/elementSubRegions/cellBlock01/elementVolume', (1,), float)])
    def test_wrapper_read(self, target_key, target_shape, target_type, geosx_problem):
        """
        Test whether wrappers can be read, and whether their shape and value types are correct

        Args:
            target_key (str): Key value for the target wrapper
            target_shape (tuple): Expected shape for the target wrapper
            target_type (dtype): Expected type for the wrapper elements
            geosx_proplem (pygeosx.problem): The problem handle
        """
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

    @pytest.mark.parametrize('target_key',
                             [('Events/time'),
                              ('Events/cycle'),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ReferencePosition'),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ghostRank'),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/ElementRegions/elementRegionsGroup/dummy/elementSubRegions/cellBlock01/elementVolume')])
    def test_wrapper_write(self, target_key, geosx_problem):
        """
        Test writing to the wrappers

        Args:
            target_key (str): Key value for the target wrapper
            target_shape (tuple): Expected shape for the target wrapper
            target_type (dtype): Expected type for the wrapper elements
            geosx_proplem (pygeosx.problem): The problem handle
        """
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

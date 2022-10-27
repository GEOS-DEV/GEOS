
import pylvarray
import pygeosx
import numpy as np
import pytest
from pygeosx_tools import wrapper


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
        rank = 0
        args = ('pygeosx', '-i', 'pygeosx_only.xml')
        problem = pygeosx.initialize(rank, args)
        pygeosx.apply_initial_conditions()
        # print(wrapper.get_matching_wrapper_path(problem, ['cellBlock01']))
        return problem

    @pytest.mark.parametrize('target_key, target_shape, target_type',
                             [('Events/time', (1,), float),
                              ('Events/cycle', (1,), np.int32),
                              ('Events/python/target', 21, str),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ReferencePosition', (8, 3), float),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/nodeManager/ghostRank', (8,), np.int32),
                              ('domain/MeshBodies/mesh1/meshLevels/Level0/ElementRegions/elementRegionsGroup/dummy/elementSubRegions/cellBlock01/elementVolume', (1,), float)])
    def test_wrapper_eval(self, target_key, target_shape, target_type, geosx_problem):
        """
        Test the shape and value of wrappers

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

        # print('test')
        # print(s, target_shape)
        # print(type(val), target_type)
        assert s == target_shape
        assert isinstance(val, target_type)



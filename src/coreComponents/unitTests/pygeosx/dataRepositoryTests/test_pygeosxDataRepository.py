
import pylvarray
import pygeosx
import numpy as np
import os
import sys
import pytest
from mpi4py import MPI
from geosx_xml_tools.main import preprocess_parallel
from pygeosx_tools import wrapper


class TestStuff(object):

    @pytest.fixture(scope='class')
    def geosx_problem(self):
        target_xml = 'sedov_pygeosx.xml'
        sys.argv = ['pygeosx', '-i', target_xml, '-c', 'sedov_pygeosx_preprocessed.xml']
        args = preprocess_parallel()

        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
        problem = pygeosx.initialize(rank, args)
        pygeosx.apply_initial_conditions()
        return problem

    def test_get_2D_array_float(self, geosx_problem):
        location_key = 'domain/MeshBodies/mesh1/meshLevels/Level0/ElementRegions/elementRegionsGroup/Region2/elementSubRegions/cb1/elementCenter'
        x = wrapper.get_wrapper(geosx_problem, location_key)
        assert isinstance(x, np.ndarray)
        assert np.shape(x) == (1000, 3)
        assert isinstance(x[0, 0], float)

    def test_get_1D_array_int(self, geosx_problem):
        ghost_key = 'domain/MeshBodies/mesh1/meshLevels/FE1/ElementRegions/elementRegionsGroup/Region2/elementSubRegions/cb1/ghostRank'
        x = wrapper.get_wrapper(geosx_problem, ghost_key)
        assert isinstance(x, np.ndarray)
        assert isinstance(x[0], np.int32)


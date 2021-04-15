
import sys
import numpy as np
from mpi4py import MPI
import pygeosx
from pygeosx_tools import wrapper


def surface_mapping(x, z_transition=0.5, z_surf=1.0):
    """
    Function to deform the mesh to match some surface topography.
    In this example, the surface elevation is set to:
    0.1 * sin(2 * pi * r / max(r)), and the mesh deformation
    is tapered as a function of depth.

    Args:
        x (np.ndarray): the object locations
        z_transition (float): depth of the zone to be deformed
        z_surf (float): average location of the model surface

    Returns:
        np.ndarray: new coordinates
    """
    R = x[:, 0]**2 + x[:, 1]**2
    delta_surf = 0.1 * np.sin(2.0 * np.pi * R / np.amax(R))
    delta_z = delta_surf * np.maximum(x[:, 2] - z_surf + z_transition, 0.0) / z_transition
    return x[:, 2] + delta_z


def run_problem():
    """
    Run the GEOSX problem
    """
    # Initialize the code
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    problem = pygeosx.initialize(rank, sys.argv)

    # Rather than specifying the wrapper paths explicitly,
    # search for them using a set of filters
    location_keys = []
    location_keys.append(wrapper.get_matching_wrapper_path(problem, ['ReferencePosition']))
    location_keys.append(wrapper.get_matching_wrapper_path(problem, ['Region2', 'elementCenter']))

    # Map the node, element locations
    print('Warning: this mapping will change the mesh, but for things to work correctly')
    print('         we would need to re-calculate the geometric parameters, shape functions')
    for ka in location_keys:
        wrapper.set_wrapper_with_function(problem, ka, ka, surface_mapping, target_index=2)

    # Apply initial conditions
    pygeosx.apply_initial_conditions()

    # Run the code
    while pygeosx.run() != pygeosx.COMPLETED:
        # Do something
        pass


if __name__ == '__main__':
    run_problem()




import numpy as np
from mpi4py import MPI
import pygeosx
from pygeosx_tools import wrapper
from geosx_xml_tools.main import preprocess_parallel


# PYGEOSX_STRESS_FN
def stress_fn(x):
    """
    Function to set stress values

    Args:
        x (np.ndarray) the element centers

    Returns:
        np.ndarray: stress values
    """
    R = x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2
    return np.sin(2.0 * np.pi * R / np.amax(R))
# PYGEOSX_STRESS_FN_END


def run_problem():
    """
    Run the GEOSX problem
    """
    # PYGEOSX_INITIALIZATION
    # Get the MPI rank
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Initialize the code and set initial conditions
    args = preprocess_parallel()
    problem = pygeosx.initialize(rank, args)
    pygeosx.apply_initial_conditions()

    # Rather than specifying the wrapper paths explicitly,
    # search for them using a set of filters
    location_key = wrapper.get_matching_wrapper_path(problem, ['Region2', 'elementCenter'])
    stress_key = wrapper.get_matching_wrapper_path(problem, ['Region2', 'shale', 'stress'])
    ghost_key = wrapper.get_matching_wrapper_path(problem, ['Region2', 'cb1', 'ghostRank'])
    # PYGEOSX_INITIALIZATION_END

    # PYGEOSX_STRESS_IC
    # Print initial stress
    wrapper.print_global_value_range(problem, stress_key, 'stress')

    # Zero out stress
    wrapper.set_wrapper_to_value(problem, stress_key, 0.0)
    wrapper.print_global_value_range(problem, stress_key, 'stress')

    # Set stress via a function
    wrapper.set_wrapper_with_function(problem, stress_key, location_key, stress_fn, target_index=0)
    wrapper.set_wrapper_with_function(problem, stress_key, location_key, stress_fn, target_index=1)
    wrapper.set_wrapper_with_function(problem, stress_key, location_key, stress_fn, target_index=2)
    wrapper.print_global_value_range(problem, stress_key, 'stress')
    # PYGEOSX_STRESS_IC_END

    # PYGEOSX_MAIN_LOOP
    # Run the code
    while pygeosx.run() != pygeosx.COMPLETED:
        wrapper.print_global_value_range(problem, stress_key, 'stress')

        # Gather/allgather tests
        tmp = wrapper.gather_wrapper(problem, stress_key)
        print(wrapper.rank, 'gather', np.shape(tmp), flush=True)

        tmp = wrapper.allgather_wrapper(problem, stress_key)
        print(wrapper.rank, 'allgather', np.shape(tmp), flush=True)

        tmp = wrapper.allgather_wrapper(problem, stress_key, ghost_key=ghost_key)
        print(wrapper.rank, 'allgather_ghost_filtered', np.shape(tmp), flush=True)
    # PYGEOSX_MAIN_LOOP_END


if __name__ == '__main__':
    run_problem()



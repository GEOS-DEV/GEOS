
import sys
import numpy as np
from mpi4py import MPI
import pygeosx
from pygeosx_tools import wrapper


def stress_fn(x):
  """
  @brief function to set stress
  @param x the element centers
  @returns a 1D vector of stress values
  """
  R = x[:, 0]**2 + x[:, 1]**2 + x[:, 2]**2
  return np.sin(2.0 * np.pi * R / np.amax(R))


def run_problem():
  """
  @brief Run the GEOSX problem
  """
  # Initialize the code
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  problem = pygeosx.initialize(rank, sys.argv)

  # Apply initial conditions
  pygeosx.apply_initial_conditions()

  # Rather than specifying the wrapper paths explicitly,
  # search for them using a set of filters
  location_key = wrapper.get_matching_wrapper_path(problem, ['Region2', 'elementCenter'])
  stress_key = wrapper.get_matching_wrapper_path(problem, ['Region2', 'shale', 'stress'])

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

  # Run the code
  while pygeosx.run() != pygeosx.COMPLETED:
    wrapper.print_global_value_range(problem, stress_key, 'stress')


if __name__ == '__main__':
  run_problem()



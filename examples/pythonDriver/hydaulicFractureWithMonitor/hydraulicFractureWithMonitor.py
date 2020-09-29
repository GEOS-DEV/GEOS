
import sys
from mpi4py import MPI
import matplotlib.pyplot as plt
import pygeosx
from pygeosx_tools import wrapper, xml


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def run_problem():
  """
  @brief Run the GEOSX problem
  """
  # Initialize the code
  xml.apply_xml_preprocessor()
  problem = pygeosx.initialize(rank, sys.argv)

  # Apply initial conditions
  pygeosx.apply_initial_conditions()

  # Rather than specifying the wrapper paths explicitly,
  # search for them using a set of filters
  fracture_location_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'elementCenter'])
  fracture_aperture_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'effectiveAperture'])

  # Note: we will be plotting our results.
  # This will modify the fonts so that they a bit easier to read
  plot_font = {'weight': 'normal', 'size': 8}
  plt.rc('font', **plot_font)

  # Setup values to record
  records = {fracture_location_key: {'label': 'Fracture Extents (m)',
                                     'scale': 1.0,
                                     'history': [],
                                     'fhandle': plt.figure()},
             fracture_aperture_key: {'label': 'Aperture (mm)',
                                     'scale': 1e3,
                                     'history': [],
                                     'fhandle': plt.figure()},
             'time': {'label': 'Time (min)',
                      'scale': 1.0 / 60.0,
                      'history': []}}

  # Run the code
  while pygeosx.run() != pygeosx.COMPLETED:
    wrapper.run_queries(problem, records)
    wrapper.plot_history(records)


if __name__ == '__main__':
  run_problem()



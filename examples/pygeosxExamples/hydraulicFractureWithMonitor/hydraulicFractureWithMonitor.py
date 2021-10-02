
# PYGEOSX_SETUP
import sys
from mpi4py import MPI
import pygeosx
from pygeosx_tools import wrapper, xml
import matplotlib.pyplot as plt
# PYGEOSX_SETUP_END


def run_problem(plot_frequency=1):
    """
    Run the GEOSX problem
    """
    # PYGEOSX_INITIALIZATION
    # Get the MPI rank
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Initialize the code and set initial conditions
    xml.apply_xml_preprocessor()
    problem = pygeosx.initialize(rank, sys.argv)
    pygeosx.apply_initial_conditions()
    # PYGEOSX_INITIALIZATION_END

    # PYGEOSX_KEY_SEARCH
    # Rather than specifying the wrapper paths explicitly,
    # search for them using a set of filters
    fracture_location_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'elementCenter'])
    fracture_aperture_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'effectiveAperture'])
    # PYGEOSX_KEY_SEARCH_END

    # Note: we will be plotting our results.
    # This will modify the fonts so that they a bit easier to read
    plot_font = {'weight': 'normal', 'size': 8}
    plt.rc('font', **plot_font)

    # PYGEOSX_QUERY_SETUP
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
    # PYGEOSX_QUERY_SETUP_END

    # PYGEOSX_MAIN_LOOP
    # Run the code
    query_count = 0
    while pygeosx.run() != pygeosx.COMPLETED:
        wrapper.run_queries(problem, records)
        query_count += 1
        if (query_count % plot_frequency == 0):
            wrapper.plot_history(records)
    # PYGEOSX_MAIN_LOOP_END


if __name__ == '__main__':
    run_problem(plot_frequency=4)



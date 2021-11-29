
import sys
from mpi4py import MPI
import pygeosx
from pygeosx_tools import wrapper
import numpy as np
import matplotlib.pyplot as plt


def run_problem():
    """
    Run the GEOSX problem
    """
    # Initialize the code
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    problem = pygeosx.initialize(rank, sys.argv)

    # Rather than specifying the wrapper paths explicitly,
    time_key = 'Events/time'
    bc_key = wrapper.get_matching_wrapper_path(problem, ['sourceTerm', 'scale'])
    fracture_location_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'elementCenter'])
    fracture_aperture_key = wrapper.get_matching_wrapper_path(problem, ['Fracture', 'effectiveAperture'])

    # Apply initial conditions
    pygeosx.apply_initial_conditions()

    # Setup monitor
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
        time = problem.get_wrapper(time_key).value(False)
        current_flow_rate = -5.0 - np.sin(np.pi * time / 60.0)
        wrapper.set_wrapper_to_value(problem, bc_key, current_flow_rate)
        if (wrapper.rank == 0):
            print('t = %1.4f, q = %1.4f' % (time, current_flow_rate), flush=True)

        wrapper.run_queries(problem, records)
        wrapper.plot_history(records, show_figures=True)


if __name__ == '__main__':
    run_problem()



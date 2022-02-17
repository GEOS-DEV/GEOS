
import sys
from mpi4py import MPI
import pygeosx
from pygeosx_tools import xml


def run_problem():
    """
    Run the GEOSX problem
    """
    # PYGEOSX_INITIALIZATION
    # Get the MPI rank
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Preprocess the xml file
    xml.apply_xml_preprocessor()
    
    # Initialize the code and set initial conditions
    problem = pygeosx.initialize(rank, sys.argv)
    pygeosx.apply_initial_conditions()
    while pygeosx.run() != pygeosx.COMPLETED:
        pass

if __name__ == '__main__':
    run_problem()



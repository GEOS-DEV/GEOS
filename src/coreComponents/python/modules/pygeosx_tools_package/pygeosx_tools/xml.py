import sys
from mpi4py import MPI
from geosx_xml_tools import xml_processor


# Get the MPI rank
comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def apply_xml_preprocessor():
    """
    Applies the xml preprocessor to the input file
    before handing it to GEOSX, and modifies the input
    arguments to point to the new file
    """
    new_input_file = ''
    file_index = sys.argv.index('-i') + 1
    if (rank == 0):
        print('Applying xml preprocessor...')
        new_input_file = xml_processor.process(sys.argv[file_index], '%s.processed' % (sys.argv[file_index]))
        print('  the compiled filename is: %s' % (new_input_file))

    # Broadcast and set the new input file name
    sys.argv[file_index] = comm.bcast(new_input_file, root=0)


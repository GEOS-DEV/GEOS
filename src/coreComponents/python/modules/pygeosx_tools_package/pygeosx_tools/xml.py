from geosx_xml_tools import xml_processor
from pygeosx_tools import parallel_io


def apply_xml_preprocessor(geosx_args):
    """
    Applies the xml preprocessor to the input file
    before handing it to GEOSX, and modifies the input
    arguments to point to the new file

    Args:
        geosx_args (list): User arguments supplied to GEOSX
    """
    new_input_file = ''
    file_index = geosx_args.index('-i') + 1
    if (parallel_io.rank == 0):
        print('Applying xml preprocessor...', flush=True)
        new_input_file = xml_processor.process(geosx_args[file_index], '%s.processed' % (geosx_args[file_index]))
        print('  the compiled filename is: %s' % (new_input_file), flush=True)

    # Broadcast and set the new input file name
    geosx_args[file_index] = parallel_io.comm.bcast(new_input_file, root=0)


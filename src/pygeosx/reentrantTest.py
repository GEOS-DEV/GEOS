import sys
import os

from mpi4py import MPI

import pygeosx


def main():
    # Get the MPI rank
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        print('In python')
        sys.stdout.flush()

    # If not a restart run then we'll run the SSLE-sedov problem first.
    if '-r' not in sys.argv:
        arg_copy = list(sys.argv)

        xml_index = arg_copy.index('-i') + 1
        arg_copy[xml_index] = os.path.join(
            os.path.dirname(__file__), '..', '..',
            'integratedTests', 'update', 'run', 'solidMechanicsSSLE', 'SSLE-sedov.xml'
        )

        output_index = arg_copy.index('-o') + 1
        arg_copy[output_index] = sys.argv[output_index] + '_ssle_dummy'

        problem = pygeosx.initialize(rank, arg_copy)
        pygeosx.apply_initial_conditions()

        # run to completion
        while pygeosx.run() != pygeosx.COMPLETED:
            pass

        if rank == 0:
            print('\n\nIn python second time around')
            sys.stdout.flush()

        problem = pygeosx.reinit(sys.argv)
    else:
        problem = pygeosx.initialize(rank, sys.argv)

    pygeosx.apply_initial_conditions()

    while pygeosx.run() != pygeosx.COMPLETED:
        pass


if __name__ == '__main__':
    main()

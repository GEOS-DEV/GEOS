import pygeosx
from mpi4py import MPI
import sys

from AcousticShot import *
from shotFileManager import *

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def main():
    if len(sys.argv) == 5:
        shot_file = sys.argv[3]
        tracePath = sys.argv[4]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:3], shot_list, tracePath)

    elif len(sys.argv) == 7:
        shot_file = sys.argv[5]
        tracePath = sys.argv[6]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:5], shot_list, tracePath)

    elif len(sys.argv) == 9:
        shot_file = sys.argv[7]
        tracePath = sys.argv[8]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:7], shot_list, tracePath)

    elif len(sys.argv) == 11:
        shot_file = sys.argv[9]
        tracePath = sys.argv[10]

        shot_list = readShotList(shot_file)

        shot_simul(rank, sys.argv[0:9], shot_list, tracePath)



def shot_simul(rank, xml, shot_list, tracePath):
    """
    Parameters
    ----------

    shot_list : list
        A list containing sets of Source and ReceiverSet objects

    dt : float
        Time step for the solver
    """

    problem = initialize(rank, xml)
    acoustic_shots(rank, problem, shot_list, tracePath)



def initialize(rank, xml):
    """ Grouping of pygeox initializations

    Return
    ------
    problem :
        The pygeosx ProblemManager

    Notes
    -----
    Need to give MPI rank at this point for initialization.
    Conflict with first initialization to get the list of shots
    """

    problem = pygeosx.initialize(rank, xml)
    pygeosx.apply_initial_conditions()

    return problem


if __name__ == "__main__":
    main()

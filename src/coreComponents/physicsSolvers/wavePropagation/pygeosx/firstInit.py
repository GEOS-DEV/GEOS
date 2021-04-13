import pygeosx
from mpi4py import MPI
import sys

from mesh import calculDt, domainBoundary
from shotFileManager import exportInitVariable

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():

    problem = pygeosx.initialize(0, sys.argv)
    pygeosx.apply_initial_conditions()
    maxT = problem.get_wrapper("Events/maxTime").value()[0]
    dt = calculDt(problem)
    boundary_box  = domainBoundary(problem)

    exportInitVariable(maxT, dt, boundary_box)


if __name__ == "__main__":
    main()

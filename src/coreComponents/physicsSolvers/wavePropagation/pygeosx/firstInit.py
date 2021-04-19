import pygeosx
from mpi4py import MPI
import sys

from mesh import calculDt, domainBoundary
from shotFileManager import exportInitVariable

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def main():

    problem = pygeosx.initialize(rank, sys.argv)
    pygeosx.apply_initial_conditions()
    
    if rank==0:
        maxT = problem.get_wrapper("Events/maxTime").value()[0]
        dt = calculDt(problem)
        boundary_box  = domainBoundary(problem)

        exportInitVariable(maxT, dt, boundary_box)


if __name__ == "__main__":
    main()

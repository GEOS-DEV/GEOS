import pygeosx
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

def test_forwardPropagationAcoustic():
    argv = ["testPygeosxAcoustic.py", "-i", "acous_5x5x5.xml"]
    solver = "acousticSolver"
    forwardPropagation(argv, solver)

def forwardPropagation(argv, solver)
    problem = pygeosx.initialize(rank, argv)
    pygeosx.apply_initial_conditions()

    solver = problem.get_group(f"/Solvers/{solver}")

    t=0
    dt = 0.005
    tmax = 0.1
    cycle = 0

    while t < tmax:
        solver.execute(t, dt)
        t+=dt
        cycle+=1

    rcvs_pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()

if __name__ == "__main__":
    test_forwardPropagationAcoustic()

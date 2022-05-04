import pygeosx
from mpi4py import MPI
import numpy as np
import sys

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

problem = pygeosx.initialize(rank, sys.argv)
pygeosx.apply_initial_conditions()

solver = problem.get_group("/Solvers/acousticSolver")

t=0
dt = 0.005
tmax = 1.5
cycle = 0

while t < tmax:
    solver.execute(t, dt)
    t+=dt
    print("time = %.3fs," % t, "dt = %.4f," % dt, "cycle =", cycle)
    cycle+=1

rcvs_pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()

print("\npressure at receivers = \n", rcvs_pressure)

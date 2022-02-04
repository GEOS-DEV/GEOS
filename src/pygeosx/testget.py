import pygeosx
import sys
import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

problem = pygeosx.initialize(rank, sys.argv)

solver = problem.get_group("/Solvers/acousticSolver")

collection = problem.get_group("/Tasks/waveFieldCollection")
output = problem.get_group("Outputs/waveFieldOutput")

pygeosx.apply_initial_conditions()
time = 0
dt = 0.005
for i in range(100):
    solver.execute(time, dt)
    if i==99 :
        collection.collect(time, dt)
        output.output(time, dt)

    time += dt
    print("time =", time)
pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()


pygeosx._finalize()

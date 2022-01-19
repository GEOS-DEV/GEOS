import pygeosx
import sys
import h5py
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
    solver.solverStep(time, dt)
    if i==99 :
        collection.collect(time, dt)
        output.output(time, dt)

    time += dt
    print("time =", time)
pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()
print(pressure)


filename = "waveField.hdf5"

if rank == 0:
    with h5py.File(filename, "r") as f:
        a_group_key = list(f.keys())[0]

        # Get the data
        data = list(f[a_group_key])

    print(data[0])
pygeosx._finalize()

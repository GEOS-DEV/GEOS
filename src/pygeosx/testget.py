import pygeosx
import sys
import h5py
import numpy as np

problem = pygeosx.initialize(0, sys.argv)


solver = problem.get_group("/Solvers/acousticSolver")

collection = problem.get_group("/Tasks/waveFieldCollection")
output = problem.get_group("Outputs/waveFieldOutput")

pygeosx.apply_initial_conditions()
for t in range(1000):
    solver.solverStep(t, 0.005)
    if t==999:
        collection.collect(t, 0.005)
        output.output(t, 0.005)

pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()
print(pressure)


filename = "waveField.hdf5"

with h5py.File(filename, "r") as f:
    a_group_key = list(f.keys())[0]

    # Get the data
    data = list(f[a_group_key])

print(data[0])

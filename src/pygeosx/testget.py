import pygeosx
import sys

problem = pygeosx.initialize(0, sys.argv)


solver = problem.get_group("/Solvers/acousticSolver")

#history = pygeosx.pyhistory.History()

pygeosx.apply_initial_conditions()
solver.solverStep(0, 0.005)
#history.collect("waveField",0,0.005)

pressure = solver.get_wrapper("pressureNp1AtReceivers").value().to_numpy()
print(pressure)

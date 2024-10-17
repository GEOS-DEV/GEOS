# ------------------------------------------------------------------------------------------------------------
# SPDX-License-Identifier: LGPL-2.1-only
#
# Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
# Copyright (c) 2018-2024 Total, S.A
# Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
# Copyright (c) 2023-2024 Chevron
# Copyright (c) 2019-     GEOS/GEOSX Contributors
# Copyright (c) 2019-     INRIA project-team Makutu 
# All rights reserved
#
# See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
# ------------------------------------------------------------------------------------------------------------

import numpy as np
from mpi4py import MPI

from utilities.input import XML
from utilities.solvers import ReservoirSolver
from darts_model import Model 

def run_darts_model(xml_name: str, darts_model=None):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    xml = XML(xml_name)

    solver = ReservoirSolver()
    solver.initialize(rank=rank, xml=xml)
    solver.initialDt = solver.geosx.get_wrapper("Events/solver_1/forceDt").value()[0]
    solver.dt = solver.initialDt

    functions = solver.geosx.get_group("/Functions").groups()
    for func in functions:
        if hasattr(func, 'setAxes') and darts_model is not None:
          func.setAxes( darts_model.physics.n_vars, 
                        darts_model.physics.n_ops, 
                        list(darts_model.physics.axes_min), 
                        list(darts_model.physics.axes_max), 
                        list(darts_model.physics.n_axes_points) )
          func.setEvaluateFunction(darts_model.physics.reservoir_operators[0].evaluate)
          print("Adaptive OBL interpolator is configured.")

    if rank == 0:
       print()
       print("===============================================================")
       print(" Running pygeosx with solver:", solver.solver)
       print(" simulation duration: ", solver.maxTime)
       print(" initial timestep: ", solver.initialDt)
       print("===============================================================")
    solver.applyInitialConditions()

    time = 0
    cycle = 0
   
    solver.outputVtk(time)
    while time < solver.maxTime:
        if rank == 0:
          if solver.initialDt is not None:
              print(f"time = {time:.3f}s, dt = {solver.initialDt:.4f}, iter = {cycle+1}")
          else:
              print(f"time = {time:.3f}s, dt is None, iter = {cycle+1}")
        solver.execute(time)
        solver.updateTimeStep()
        solver.outputVtk(time)
        time += solver.initialDt
        cycle += 1
    solver.cleanup(time)
    
    comm.Barrier()


# run adaptive OBL
print("\n" + "="*30 + " RUNNING ADAPTIVE OBL " + "="*30 + "\n")
darts_model = Model()
run_darts_model(xml_name="input_file_adaptive.xml", darts_model=darts_model)

# run static OBL
print("\n" + "="*30 + " RUNNING STATIC OBL " + "="*30 + "\n")
run_darts_model(xml_name="input_file_static.xml")

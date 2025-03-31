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

from geos.pygeos_tools.utilities.input import XML
from geos.pygeos_tools.utilities.solvers import ReservoirSolver

from model import Model

def run_darts_model(domain: str, xml_name: str, darts_model=None):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    xml = XML(xml_name)

    solver = ReservoirSolver()
    solver.initialize(rank=rank, xml=xml)

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

    solver.applyInitialConditions()

    time = 0
    cycle = 0
    solver.initialDt = 8.64 
    solver.dt = solver.initialDt
    
    solver.outputVtk(time)
    while time < solver.maxTime:
        # choose new timestep
        if domain == '1D':
            if time < 48: solver.initialDt = 4.0
            elif time < 240: solver.initialDt = 8.64
            elif time < 3600: solver.initialDt = 86.4
            elif time < 6 * 8640: solver.initialDt = 240.0
            elif time < 2 * 86400: solver.initialDt = 900.0
            else: solver.initialDt = 3600.0
        elif domain == '2D':
            if time < 24: solver.initialDt = 4.0
            elif time < 120: solver.initialDt = 8.64
            elif time < 300: solver.initialDt = 2 * 8.64
            elif time < 1400: solver.initialDt = 60.0
            elif time < 4 * 3600: solver.initialDt = 300.0
            elif time < 9 * 3600: solver.initialDt = 200.0
            else: solver.initialDt = 100.0
        if rank == 0:
            print(f"time = {time:.3f}s, dt = {solver.initialDt:.4f}, step = {cycle+1}")
        # run simulation
        solver.execute(time)
        time += solver.initialDt
        if cycle % 5 == 0:
            solver.outputVtk(time)
        cycle += 1
    solver.cleanup(time)
    
    comm.Barrier()


darts_model = Model()
run_darts_model(domain='2D', xml_name="2d_setup.xml", darts_model=darts_model)
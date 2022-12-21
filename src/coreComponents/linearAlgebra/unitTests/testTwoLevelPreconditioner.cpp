/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testKrylovSolvers.cpp
 */

#include "common/DataTypes.hpp"
#include "linearAlgebra/solvers/PreconditionerTwoLevel.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"

#include <gtest/gtest.h>

using namespace geosx;

int main( int argc, char * * argv )
{
  geosx::testing::LinearAlgebraTestScope scope( argc, argv );
  
  HypreMatrix Ah;
  HypreMatrix AH;

  Ah.read("matrixFine", MPI_COMM_GEOSX );
  AH.read("matrixCoarse", Ah.comm() );

  GEOSX_LOG_RANK_VAR( Ah.numLocalRows() );
  GEOSX_LOG_RANK_VAR( AH.numLocalRows() );

  PreconditionerTwoLevel< HypreInterface > precond( Ah, AH);

  // Create random tue solution vector
  HypreVector sol_true;
  sol_true.create( Ah.numLocalCols(), Ah.comm() );
  sol_true.rand( 1984 );
GEOSX_LOG_RANK("true sol DONE");

  // Create and compute the rhs vector
  HypreVector rhs;
  rhs.create( Ah.numLocalRows(), Ah.comm() );
  Ah.apply( sol_true, rhs );
GEOSX_LOG_RANK("rhs computed");

  // Create and zero out the computed solution vector
  HypreVector sol_comp;
  sol_comp.create( sol_true.localSize(), sol_true.comm() );
  sol_comp.zero();
GEOSX_LOG_RANK("iterative sol INITIALIZED");

  // Create PCG solver
  LinearSolverParameters params;
  params.isSymmetric = true;
  params.logLevel = 2;
  params.solverType = LinearSolverParameters::SolverType::cg;
  std::unique_ptr< KrylovSolver< HypreVector > > const solver = KrylovSolver< HypreVector >::create( params, Ah, precond );

  // Solve linear system
  solver->solve( rhs, sol_comp );

  // Check solution
  EXPECT_TRUE( solver->result().success() );

  HypreVector sol_diff( sol_comp );
  sol_diff.axpy( -1.0, sol_true );
  GEOSX_LOG_RANK_VAR( params.krylov.relTolerance );
  GEOSX_LOG_RANK_VAR( sol_diff.norm2() / sol_true.norm2() );

  GEOSX_LOG_RANK_VAR("Completed");
}

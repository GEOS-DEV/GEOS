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
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerTwoLevel.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "linearAlgebra/solvers/CgSolver.hpp"
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

  // Create unit rhs
  HypreVector rhs;
  rhs.create( Ah.numLocalRows(), Ah.comm() );
  rhs.set(1.0);

  // Create solution vector initialized to zero
  HypreVector sol_comp;
  sol_comp.create( rhs.localSize(), rhs.comm() );
  sol_comp.zero();

  // Create preconditioner and PCG solver
  LinearSolverParameters params;
  params.isSymmetric = true;
  params.logLevel = 2;
  params.solverType = LinearSolverParameters::SolverType::cg;
  
  PreconditionerTwoLevel< HypreInterface > precond( Ah, AH);
  //PreconditionerJacobi< HypreInterface > precond;
  //precond.setup(Ah);

  CgSolver< HypreVector > solver( params, Ah, precond );

  // Solve linear system
  solver.solve( rhs, sol_comp );

  // Check solution
  EXPECT_TRUE( solver.result().success() );
}

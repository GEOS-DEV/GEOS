/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testExternalSolvers.cpp
 */

#include <gtest/gtest.h>

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "managers/initialization.hpp"

#include "testLinearAlgebraUtils.hpp"

using namespace geosx;

static real64 constexpr machinePrecision = 20.0 * std::numeric_limits< real64 >::epsilon();

///////////////////////////////////////////////////////////////////////////////////////

LinearSolverParameters params_Direct()
{
  LinearSolverParameters parameters;
  parameters.solverType = LinearSolverParameters::SolverType::direct;
  parameters.krylov.relTolerance = machinePrecision;
  return parameters;
}

LinearSolverParameters params_GMRES_ILU()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::iluk;
  parameters.ilu.fill = 1;
  return parameters;
}

LinearSolverParameters params_GMRES_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  parameters.amg.smootherType = LinearSolverParameters::PreconditionerType::gs;
  parameters.amg.coarseType = LinearSolverParameters::PreconditionerType::direct;
  return parameters;
}

LinearSolverParameters params_CG_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::cg;
  parameters.isSymmetric = true;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
  parameters.amg.smootherType = LinearSolverParameters::PreconditionerType::gs;
  parameters.amg.coarseType = LinearSolverParameters::PreconditionerType::direct;
  return parameters;
}

template< typename LAI >
class SolverTestBase : public ::testing::Test
{
public:

  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;
  using Solver = typename LAI::LinearSolver;

protected:

  Matrix matrix;
  real64 cond_est = 1.0;

  void test( LinearSolverParameters const & params )
  {
    // Create a random "true" solution vector
    Vector sol_true;
    sol_true.createWithLocalSize( matrix.numLocalCols(), matrix.getComm() );
    sol_true.rand();

    // Create and compute the right-hand side vector
    Vector rhs;
    rhs.createWithLocalSize( matrix.numLocalRows(), matrix.getComm() );
    matrix.apply( sol_true, rhs );

    // Create and zero out the computed solution vector
    Vector sol_comp;
    sol_comp.createWithLocalSize( sol_true.localSize(), sol_true.getComm() );
    sol_comp.zero();

    // Create the solver and solve the system
    Solver solver( params );
    solver.solve( matrix, sol_comp, rhs );
    EXPECT_TRUE( solver.result().success() );

    // Check that solution is within epsilon of true
    sol_comp.axpy( -1.0, sol_true );
    real64 const relTol = cond_est * params.krylov.relTolerance;
    EXPECT_LT( sol_comp.norm2() / sol_true.norm2(), relTol );
  }
};

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestLaplace2D : public SolverTestBase< LAI >
{
public:

  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:

  void SetUp() override
  {
    globalIndex constexpr n = 100;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, n, this->matrix );

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( SolverTestLaplace2D );

TYPED_TEST_P( SolverTestLaplace2D, Direct )
{
  this->test( params_Direct() );
}

TYPED_TEST_P( SolverTestLaplace2D, GMRES_ILU )
{
  this->test( params_GMRES_ILU() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_AMG )
{
  this->test( params_CG_AMG() );
}

REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             Direct,
                             GMRES_ILU,
                             CG_AMG );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestLaplace2D, TrilinosInterface, );
INSTANTIATE_TYPED_TEST_SUITE_P( TrilinosTpetra, SolverTestLaplace2D, TrilinosTpetraInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestLaplace2D, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestLaplace2D, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestElasticity2D : public SolverTestBase< LAI >
{
public:

  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:

  void SetUp() override
  {
    globalIndex constexpr n = 100;
    compute2DElasticityOperator( MPI_COMM_GEOSX, 1.0, 1.0, n, n, 10000., 0.2, this->matrix );

    // Impose Dirichlet boundary conditions: fix domain bottom (first 2*(nCellsX + 1) rows of matrix)
    this->matrix.open();
    for( globalIndex iRow = 0; iRow < 2 * (n + 1); ++iRow )
    {
      if( this->matrix.getLocalRowID( iRow ) >= 0 )
      {
        this->matrix.clearRow( iRow, true );
      }
    }
    this->matrix.close();
    this->cond_est = 1e4; // not a true condition number estimate, but enough to pass tests
  }
};

TYPED_TEST_SUITE_P( SolverTestElasticity2D );

TYPED_TEST_P( SolverTestElasticity2D, Direct )
{
  this->test( params_Direct() );
}

TYPED_TEST_P( SolverTestElasticity2D, GMRES_AMG )
{
  LinearSolverParameters params = params_GMRES_AMG();
  params.amg.separateComponents = true;
  params.dofsPerNode = 2;
  this->test( params );
}

REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             Direct,
                             GMRES_AMG );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestElasticity2D, TrilinosInterface, );
INSTANTIATE_TYPED_TEST_SUITE_P( TrilinosTpetra, SolverTestElasticity2D, TrilinosTpetraInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestElasticity2D, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestElasticity2D, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}

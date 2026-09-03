/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2016-2024 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2024 TotalEnergies
 * Copyright (c) 2018-2024 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2023-2024 Chevron
 * Copyright (c) 2019-     GEOS/GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testExternalSolvers.cpp
 */

#include "linearAlgebra/unitTests/testLinearAlgebraUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#if defined(GEOS_USE_HYPRE) && !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#endif
#if defined(GEOS_USE_HIP) && defined(GEOS_USE_HYPREDRV)
#include "linearAlgebra/interfaces/hypre/hypredrive.hpp"
#endif

#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <string>

#if defined(GEOS_USE_HIP) && defined(GEOS_USE_HYPREDRV)
#include <unistd.h>
#endif

using namespace geos;

///////////////////////////////////////////////////////////////////////////////////////

LinearSolverParameters params_DirectSerial()
{
  LinearSolverParameters parameters;
  parameters.solverType = geos::LinearSolverParameters::SolverType::direct;
  parameters.direct.parallel = 0;
  return parameters;
}

LinearSolverParameters params_DirectParallel()
{
  LinearSolverParameters parameters;
  parameters.solverType = geos::LinearSolverParameters::SolverType::direct;
  parameters.direct.parallel = 1;
  return parameters;
}

LinearSolverParameters params_GMRES_ILU()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::iluk;
#if defined(GEOS_USE_HIP)
  // HYPRE's HIP implementation supports ILU(0) on device memory; level-1
  // ILU uses a path that returns an error on the current ROCm stack.
  parameters.ifact.fill = 0;
#else
  parameters.ifact.fill = 1;
#endif
  return parameters;
}

LinearSolverParameters params_CG_SGS()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.isSymmetric = true;
  parameters.solverType = LinearSolverParameters::SolverType::cg;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::sgs;
#if defined(GEOS_USE_HIP)
  // HYPRE's HIP SGS implementation enters rocSPARSE csrsv, which is not
  // available on the ROCm stack used by the sanitizer test device. Keep the
  // test at the same CG coverage while using the device-supported relaxation.
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::l1jacobi;
#endif
  return parameters;
}

LinearSolverParameters params_GMRES_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = LinearSolverParameters::SolverType::gmres;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
#if defined(GEOS_USE_HIP)
  // HYPRE's PMIS GPU path invokes a rocPRIM radix sort that returns
  // hipErrorIllegalState on gfx10/gfx11 with the ROCm stack under test.
  // HMIS uses the device-supported hybrid coarsening path without that sort.
  parameters.amg.coarseningType = LinearSolverParameters::AMG::CoarseningType::HMIS;
#endif
#if defined(GEOS_USE_HIP)
  // Keep the AMG smoother on HYPRE's device-supported L1-Jacobi path; the
  // hybrid Gauss-Seidel variants call rocSPARSE triangular solves on HIP.
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::l1jacobi;
#else
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::fgs;
#endif
  parameters.amg.coarseType = geos::LinearSolverParameters::AMG::CoarseType::direct;
  return parameters;
}

LinearSolverParameters params_CG_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.isSymmetric = true;
  parameters.solverType = LinearSolverParameters::SolverType::cg;
  parameters.preconditionerType = LinearSolverParameters::PreconditionerType::amg;
#if defined(GEOS_USE_HIP)
  // HYPRE's PMIS GPU path invokes a rocPRIM radix sort that returns
  // hipErrorIllegalState on gfx10/gfx11 with the ROCm stack under test.
  // HMIS uses the device-supported hybrid coarsening path without that sort.
  parameters.amg.coarseningType = LinearSolverParameters::AMG::CoarseningType::HMIS;
#endif
#if defined(GEOS_USE_HIP)
  // The standard SGS implementation enters rocSPARSE csrsv on HIP. L1-Jacobi
  // exercises the same CG+AMG path without that unsupported triangular solve.
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::l1jacobi;
#else
  parameters.amg.smootherType = geos::LinearSolverParameters::AMG::SmootherType::sgs;
#endif
  parameters.amg.coarseType = geos::LinearSolverParameters::AMG::CoarseType::direct;
  return parameters;
}

#if defined(GEOS_USE_HIP) && defined(GEOS_USE_HYPREDRV)
/**
 * @brief Apply the HIP-only ILU workaround to this unit test's YAML.
 *
 * The production-generated YAML intentionally retains HYPRE's direct
 * triangular solve default. This test uses the test HYPRE build on gfx1100,
 * where that rocSPARSE analysis is unavailable, so make the exception local
 * to the test configuration.
 */
class ScopedHipIluTestConfiguration final
{
public:

  explicit ScopedHipIluTestConfiguration( LinearSolverParameters & params )
    : m_params( params )
  {}

  ~ScopedHipIluTestConfiguration()
  {
    if( !m_path.empty() )
    {
      std::remove( m_path.c_str() );
    }
  }

  bool apply()
  {
    if( m_params.solverType != LinearSolverParameters::SolverType::gmres )
    {
      return true;
    }

    if( m_params.preconditionerType != LinearSolverParameters::PreconditionerType::iluk &&
        m_params.preconditionerType != LinearSolverParameters::PreconditionerType::ilut )
    {
      return true;
    }

    hypre::hypredrive::InputArgsParseTarget target;
    if( !hypre::hypredrive::buildInputArgsParseTarget( m_params, target ) )
    {
      return false;
    }

    std::string const reorderingLine = "    reordering: 1\n";
    std::string::size_type const reordering = target.argument.find( reorderingLine );
    if( reordering == std::string::npos )
    {
      return false;
    }
    target.argument.replace( reordering,
                             reorderingLine.size(),
                             "    reordering: 0\n    tri_solve: 0\n" );

    m_path = "/tmp/geos-test-external-solvers-hip-ilu-" +
             std::to_string( static_cast< long long >( getpid() ) ) +
             ".yml";
    std::ofstream output( m_path );
    if( !output.good() )
    {
      m_path.clear();
      return false;
    }
    output << target.argument;
    if( !output.good() )
    {
      m_path.clear();
      return false;
    }

    m_params.hypredriveInputFile = Path( m_path.c_str() );
    return true;
  }

private:

  LinearSolverParameters & m_params;
  std::string m_path;
};
#endif

#if defined(GEOS_USE_HYPRE) && !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
TEST( HypreSolver, KeepsSetupDummyUntagged )
{
  struct KrylovDofLabelsGuard
  {
    ~KrylovDofLabelsGuard()
    {
      hypre::testing::clearKrylovDofLabels();
    }
  } clearKrylovDofLabels;

  HypreMatrix matrix;
  geos::testing::compute2DLaplaceOperator( MPI_COMM_GEOS, 10, matrix );

  array1d< int > labels( matrix.numLocalRows() );
  for( localIndex i = 0; i < labels.size(); ++i )
  {
    labels[i] = LvArray::integerConversion< int >( i % 2 );
  }
  hypre::testing::setKrylovDofLabels( labels.toViewConst() );

  HypreVector rhs;
  HypreVector sol;
  rhs.create( matrix.numLocalRows(), MPI_COMM_GEOS );
  rhs.set( 1.0 );
  sol.create( matrix.numLocalCols(), MPI_COMM_GEOS );
  sol.zero();

  LinearSolverParameters params;
  params.solverType = LinearSolverParameters::SolverType::fgmres;
  params.preconditionerType = LinearSolverParameters::PreconditionerType::iluk;
  params.krylov.relTolerance = 1e-8;
  params.krylov.maxIterations = 100;
  params.logLevel = 3;

  HypreSolver solver( params );
  ::testing::internal::CaptureStdout();
  solver.setup( matrix );
  solver.solve( rhs, sol );
  std::string const output = ::testing::internal::GetCapturedStdout();
  solver.clear();

  EXPECT_TRUE( solver.result().success() );
  EXPECT_EQ( output.find( "L2 norm of b0" ), std::string::npos );
}
#endif

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestBase : public ::testing::Test
{
public:

  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

protected:

  Matrix matrix;
  real64 cond_est = 1.0;

  void test( LinearSolverParameters const & params )
  {
    LinearSolverParameters solverParams = params;
#if defined(GEOS_USE_HIP) && defined(GEOS_USE_HYPREDRV)
    ScopedHipIluTestConfiguration hipIluConfiguration( solverParams );
    ASSERT_TRUE( hipIluConfiguration.apply() );
#endif

    // Create a random "true" solution vector
    Vector sol_true;
    sol_true.create( matrix.numLocalCols(), matrix.comm() );
    sol_true.rand( 1984 );

    // Create and compute the right-hand side vector
    Vector rhs;
    rhs.create( matrix.numLocalRows(), matrix.comm() );
    matrix.apply( sol_true, rhs );

    // Create and zero out the computed solution vector
    Vector sol_comp;
    sol_comp.create( sol_true.localSize(), sol_true.comm() );
    sol_comp.zero();

    // Create the solver and solve the system
    auto solver = LAI::createSolver( solverParams );
    solver->setup( matrix );
    solver->solve( rhs, sol_comp );
    EXPECT_TRUE( solver->result().success() );

    // Check that solution is within epsilon of true
    Vector sol_diff( sol_comp );
    sol_diff.axpy( -1.0, sol_true );
    real64 const relTol = cond_est * params.krylov.relTolerance;
    EXPECT_LT( sol_diff.norm2() / sol_true.norm2(), relTol );
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
    geos::testing::compute2DLaplaceOperator( MPI_COMM_GEOS, n, this->matrix );

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( SolverTestLaplace2D );

#if defined(GEOS_USE_SUITESPARSE)
TYPED_TEST_P( SolverTestLaplace2D, DirectSerial )
{
  LinearSolverParameters params = params_DirectSerial();
  params.isSymmetric = true;
  this->test( params );
}
#endif

#if !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
TYPED_TEST_P( SolverTestLaplace2D, DirectParallel )
{
  this->test( params_DirectParallel() );
}
#endif

TYPED_TEST_P( SolverTestLaplace2D, GMRES_ILU )
{
  this->test( params_GMRES_ILU() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_SGS )
{
  this->test( params_CG_SGS() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_AMG )
{
  this->test( params_CG_AMG() );
}

#if defined(GEOS_USE_CUDA) || defined(GEOS_USE_HIP)
#if defined(GEOS_USE_SUITESPARSE)
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             DirectSerial,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#endif
#else
#if defined(GEOS_USE_SUITESPARSE)
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             DirectSerial,
                             DirectParallel,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D,
                             DirectParallel,
                             GMRES_ILU,
                             CG_SGS,
                             CG_AMG );
#endif
#endif

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestLaplace2D, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestLaplace2D, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
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
    geos::testing::compute2DElasticityOperator( MPI_COMM_GEOS, 1.0, 1.0, n, n, 10000., 0.2, this->matrix );
    this->cond_est = 1e4; // not a true condition number estimate, but enough to pass tests
  }
};

TYPED_TEST_SUITE_P( SolverTestElasticity2D );

#if defined(GEOS_USE_SUITESPARSE)
TYPED_TEST_P( SolverTestElasticity2D, DirectSerial )
{
  this->test( params_DirectSerial() );
}
#endif

#if !defined(GEOS_USE_CUDA) && !defined(GEOS_USE_HIP)
TYPED_TEST_P( SolverTestElasticity2D, DirectParallel )
{
  this->test( params_DirectParallel() );
}
#endif

TYPED_TEST_P( SolverTestElasticity2D, GMRES_AMG )
{
  LinearSolverParameters params = params_GMRES_AMG();
  params.amg.separateComponents = true;
  params.dofsPerNode = 2;
  this->test( params );
}

#if defined(GEOS_USE_CUDA) || defined(GEOS_USE_HIP)
#if defined(GEOS_USE_SUITESPARSE)
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             DirectSerial,
                             GMRES_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             GMRES_AMG );
#endif
#else
#if defined(GEOS_USE_SUITESPARSE)
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             DirectSerial,
                             DirectParallel,
                             GMRES_AMG );
#else
REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D,
                             DirectParallel,
                             GMRES_AMG );
#endif
#endif

#ifdef GEOS_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestElasticity2D, TrilinosInterface, );
#endif

#ifdef GEOS_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestElasticity2D, HypreInterface, );
#endif

#ifdef GEOS_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestElasticity2D, PetscInterface, );
#endif

int main( int argc, char * * argv )
{
  geos::testing::LinearAlgebraTestScope scope( argc, argv );
  return RUN_ALL_TESTS();
}

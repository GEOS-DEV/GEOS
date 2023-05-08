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
 * @file RichardsonSolver.cpp
 */

#include "RichardsonSolver.hpp"

#include "common/TimingMacros.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geos
{

template< typename VECTOR >
RichardsonSolver< VECTOR >::RichardsonSolver( LinearSolverParameters params,
                                              LinearOperator< Vector > const & A,
                                              LinearOperator< Vector > const & M )
  : KrylovSolver< VECTOR >( std::move( params ), A, M ),
  m_omega( m_params.relaxation.weight )
{}

template< typename VECTOR >
void RichardsonSolver< VECTOR >::solve( Vector const & b,
                                        Vector & x ) const
{
  GEOS_MARK_FUNCTION;
  Stopwatch watch;

  VectorTemp r = createTempVector( b );
  VectorTemp z = createTempVector( b );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Compute the target absolute tolerance
  real64 const rnorm0 = r.norm2();
  real64 const absTol = rnorm0 * m_params.krylov.relTolerance;

  // Initialize iteration state
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.clear();

  integer & k = m_result.numIterations;
  for( k = 0; k <= m_params.krylov.maxIterations; ++k )
  {
    real64 const rnorm = r.norm2();
    m_residualNorms.emplace_back( rnorm );
    logProgress();

    // Convergence check on ||rk||/||r0||
    if( rnorm <= absTol )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // Update: z = Mr, x = x + w*z, r = b - Ax
    m_precond.apply( r, z );
    x.axpy( m_omega, z );
    m_operator.residual( x, b, r );
  }

  m_result.residualReduction = rnorm0 > 0.0 ? m_residualNorms.back() / rnorm0 : 0.0;
  m_result.solveTime = watch.elapsedTime();
  logResult();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class RichardsonSolver< TrilinosInterface::ParallelVector >;
template class RichardsonSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class RichardsonSolver< HypreInterface::ParallelVector >;
template class RichardsonSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class RichardsonSolver< PetscInterface::ParallelVector >;
template class RichardsonSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

} // geosx

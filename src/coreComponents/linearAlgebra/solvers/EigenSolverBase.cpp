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
 * @file EigenSolverBase.cpp
 */

#include "EigenSolverBase.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/solvers/ArnoldiEigenSolver.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace geos
{

template< typename VECTOR >
std::unique_ptr< GeneralizedEigenSolver< VECTOR > >
GeneralizedEigenSolver< VECTOR >::create( EigenSolverParameters const & parameters )
{
  switch( parameters.solverType )
  {
    case EigenSolverParameters::SolverType::arnoldi:
    {
      return std::make_unique< ArnoldiEigenSolver< VECTOR > >( parameters );
    }
  }
  GEOS_ERROR( "Unsupported eigensolver type" );
  return nullptr;
}

template< typename VECTOR >
VECTOR GeneralizedEigenSolver< VECTOR >::makeVector( VECTOR const & prototype )
{
  VECTOR v;
  v.create( prototype.localSize(), prototype.comm() );
  v.zero();
  return v;
}

template< typename VECTOR >
void GeneralizedEigenSolver< VECTOR >::finalizeResult( Problem const & problem,
                                                       std::vector< VECTOR > & modes,
                                                       EigenSolverResult & result ) const
{
  size_t const n = modes.size();
  result.eigenvalues.resize( n );
  result.residuals.resize( n );
  if( n == 0 )
  {
    return;
  }

  VECTOR Mx = makeVector( modes[0] );
  VECTOR Kx = makeVector( modes[0] );

  for( size_t i = 0; i < n; ++i )
  {
    // M-normalize, then evaluate the Rayleigh quotient and the residual of the original pencil
    problem.mass.apply( modes[i], Mx );
    real64 const xMx = modes[i].dot( Mx );
    GEOS_ERROR_IF( !( xMx > 0.0 ), GEOS_FMT( "Eigenvector {} has a non-positive M-norm", i ) );
    real64 const scaling = 1.0 / std::sqrt( xMx );
    modes[i].scale( scaling );
    Mx.scale( scaling );

    problem.stiffness.apply( modes[i], Kx );
    real64 const lambda = modes[i].dot( Kx );
    Kx.axpy( -lambda, Mx );

    real64 const denominator = std::fabs( lambda - m_params.shift ) * Mx.norm2();
    result.eigenvalues[i] = lambda;
    result.residuals[i] = denominator > 0.0 ? Kx.norm2() / denominator : Kx.norm2();
  }

  // Sort by ascending eigenvalue
  std::vector< size_t > order( n );
  std::iota( order.begin(), order.end(), 0 );
  std::stable_sort( order.begin(), order.end(), [&]( size_t const a, size_t const b )
  {
    return result.eigenvalues[a] < result.eigenvalues[b];
  } );

  std::vector< VECTOR > sorted;
  sorted.reserve( n );
  array1d< real64 > sortedValues( n );
  array1d< real64 > sortedResiduals( n );
  for( size_t i = 0; i < n; ++i )
  {
    sorted.push_back( std::move( modes[order[i]] ) );
    sortedValues[i] = result.eigenvalues[order[i]];
    sortedResiduals[i] = result.residuals[order[i]];
  }
  modes = std::move( sorted );
  result.eigenvalues = std::move( sortedValues );
  result.residuals = std::move( sortedResiduals );
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class GeneralizedEigenSolver< TrilinosInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_HYPRE
template class GeneralizedEigenSolver< HypreInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_PETSC
template class GeneralizedEigenSolver< PetscInterface::ParallelVector >;
#endif

} // namespace geos

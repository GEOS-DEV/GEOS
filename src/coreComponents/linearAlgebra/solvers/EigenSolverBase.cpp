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
#include "linearAlgebra/solvers/LobpcgEigenSolver.hpp"

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
    case EigenSolverParameters::SolverType::lobpcg:
    {
      return std::make_unique< LobpcgEigenSolver< VECTOR > >( parameters );
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
typename GeneralizedEigenSolver< VECTOR >::ConstraintSpace
GeneralizedEigenSolver< VECTOR >::makeConstraintSpace( Problem const & problem, VECTOR const & prototype )
{
  ConstraintSpace space;
  space.vectors.reserve( problem.constraints.size() );
  space.massVectors.reserve( problem.constraints.size() );
  for( VECTOR const * constraint : problem.constraints )
  {
    VECTOR v = makeVector( prototype );
    v.copy( *constraint );
    VECTOR mv = makeVector( prototype );
    problem.mass.apply( v, mv );
    real64 const norm0 = std::sqrt( std::max( v.dot( mv ), 0.0 ) );
    if( norm0 <= 0.0 )
    {
      continue;
    }
    space.project( v );
    problem.mass.apply( v, mv );
    real64 const norm = std::sqrt( std::max( v.dot( mv ), 0.0 ) );
    if( norm <= 1.0e-8 * norm0 )
    {
      continue;
    }
    v.scale( 1.0 / norm );
    mv.scale( 1.0 / norm );
    space.vectors.push_back( std::move( v ) );
    space.massVectors.push_back( std::move( mv ) );
  }
  return space;
}

template< typename VECTOR >
EigenSolverResult GeneralizedEigenSolver< VECTOR >::returnConstraintsOnly( Problem const & problem,
                                                                           ConstraintSpace & constraints,
                                                                           stdVector< VECTOR > & modes ) const
{
  EigenSolverResult result;
  modes.clear();
  for( integer i = 0; i < m_params.numEigenvalues && i < constraints.size(); ++i )
  {
    modes.push_back( std::move( constraints.vectors[i] ) );
  }
  result.converged = true;
  result.numConverged = LvArray::integerConversion< integer >( modes.size() );
  finalizeResult( problem, modes, result );
  return result;
}

template< typename VECTOR >
void GeneralizedEigenSolver< VECTOR >::finalizeResult( Problem const & problem,
                                                       stdVector< VECTOR > & modes,
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

  // Eigenvalue error bound ||K x - lambda M x|| / ||M x|| of each M-normalized pair
  stdVector< real64 > errorBound( n );
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

    result.eigenvalues[i] = lambda;
    errorBound[i] = Kx.norm2() / Mx.norm2();
  }

  // The bounds are relative to the spectral scale of the computed pairs and of the shift, not to |lambda - shift|:
  // the latter vanishes for the modes closest to the shift, which are the modes shift-and-invert is meant to
  // find (and the rigid-body modes for a zero shift), and would turn their rounding error into a residual of order one.
  real64 spectralScale = std::fabs( m_params.shift );
  for( size_t i = 0; i < n; ++i )
  {
    spectralScale = std::max( spectralScale, std::fabs( result.eigenvalues[i] ) );
  }
  for( size_t i = 0; i < n; ++i )
  {
    result.residuals[i] = spectralScale > 0.0 ? errorBound[i] / spectralScale : errorBound[i];
  }

  // Sort by ascending eigenvalue
  stdVector< size_t > order( n );
  std::iota( order.begin(), order.end(), 0 );
  std::stable_sort( order.begin(), order.end(), [&]( size_t const a, size_t const b )
  {
    return result.eigenvalues[a] < result.eigenvalues[b];
  } );

  stdVector< VECTOR > sorted;
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

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
 * @file LobpcgEigenSolver.cpp
 */

#include "LobpcgEigenSolver.hpp"

#include "common/Stopwatch.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

#include <algorithm>
#include <cmath>

namespace geos
{

namespace
{

using DenseMatrix = array2d< real64, MatrixLayout::COL_MAJOR_PERM >;

/// Relative threshold on the eigenvalues of the scaled Gram matrix below which a direction is dropped
real64 constexpr svqbTolerance = 1.0e-8;

/**
 * @brief Rayleigh-Ritz extraction on the span of the vectors S for the pencil (K, M).
 * @tparam VECTOR type of the vectors
 * @param[in] S basis vectors
 * @param[in] KS images of the basis vectors by K
 * @param[in] MS images of the basis vectors by M
 * @param[in] wanted number of lowest Ritz pairs requested
 * @param[out] coefficients (size(S) x wanted) coefficients of the Ritz vectors on S, M-orthonormal
 * @param[out] theta the wanted Ritz values, ascending
 *
 * The basis is M-orthonormalized with SVQB: the Gram matrix is scaled to a unit diagonal, diagonalized, and the
 * directions with a negligible eigenvalue are dropped. The Ritz problem is then a standard symmetric eigenproblem.
 */
template< typename VECTOR >
void rayleighRitz( std::vector< VECTOR * > const & S,
                   std::vector< VECTOR * > const & KS,
                   std::vector< VECTOR * > const & MS,
                   integer const wanted,
                   DenseMatrix & coefficients,
                   std::vector< real64 > & theta )
{
  integer const m = LvArray::integerConversion< integer >( S.size() );

  DenseMatrix G( m, m );
  DenseMatrix A( m, m );
  for( integer j = 0; j < m; ++j )
  {
    for( integer i = 0; i <= j; ++i )
    {
      G( i, j ) = S[i]->dot( *MS[j] );
      G( j, i ) = G( i, j );
      A( i, j ) = S[i]->dot( *KS[j] );
      A( j, i ) = A( i, j );
    }
  }

  // Scale to a unit diagonal and diagonalize
  std::vector< real64 > d( m );
  for( integer i = 0; i < m; ++i )
  {
    d[i] = G( i, i ) > 0.0 ? 1.0 / std::sqrt( G( i, i ) ) : 0.0;
  }
  DenseMatrix Gs( m, m );
  for( integer j = 0; j < m; ++j )
  {
    for( integer i = 0; i < m; ++i )
    {
      Gs( i, j ) = d[i] * G( i, j ) * d[j];
    }
  }
  array1d< real64 > sigma( m );
  DenseMatrix U( m, m );
  BlasLapackLA::matrixSymmetricEigen( Gs.toSliceConst(), sigma.toSlice(), U.toSlice() );

  std::vector< integer > kept;
  for( integer j = 0; j < m; ++j )
  {
    if( sigma[j] > svqbTolerance * sigma[m - 1] )
    {
      kept.push_back( j );
    }
  }
  integer const r = LvArray::integerConversion< integer >( kept.size() );
  GEOS_ERROR_IF( r < wanted, GEOS_FMT( "LOBPCG: the search space collapsed to {} directions, {} are needed", r, wanted ) );

  // C = D U Sigma^{-1/2} (m x r) satisfies C^T G C = I
  DenseMatrix C( m, r );
  for( integer t = 0; t < r; ++t )
  {
    real64 const s = 1.0 / std::sqrt( sigma[kept[t]] );
    for( integer i = 0; i < m; ++i )
    {
      C( i, t ) = d[i] * U( i, kept[t] ) * s;
    }
  }

  // Projected stiffness C^T A C and its eigendecomposition
  DenseMatrix AC( m, r );
  for( integer t = 0; t < r; ++t )
  {
    for( integer i = 0; i < m; ++i )
    {
      real64 sum = 0.0;
      for( integer k = 0; k < m; ++k )
      {
        sum += A( i, k ) * C( k, t );
      }
      AC( i, t ) = sum;
    }
  }
  DenseMatrix Ar( r, r );
  for( integer t = 0; t < r; ++t )
  {
    for( integer s = 0; s <= t; ++s )
    {
      real64 sum = 0.0;
      for( integer k = 0; k < m; ++k )
      {
        sum += C( k, s ) * AC( k, t );
      }
      Ar( s, t ) = sum;
      Ar( t, s ) = sum;
    }
  }
  array1d< real64 > lambda( r );
  DenseMatrix Y( r, r );
  BlasLapackLA::matrixSymmetricEigen( Ar.toSliceConst(), lambda.toSlice(), Y.toSlice() );

  coefficients.resize( m, wanted );
  theta.assign( static_cast< size_t >( wanted ), 0.0 );
  for( integer c = 0; c < wanted; ++c )
  {
    theta[c] = lambda[c];
    for( integer i = 0; i < m; ++i )
    {
      real64 sum = 0.0;
      for( integer t = 0; t < r; ++t )
      {
        sum += C( i, t ) * Y( t, c );
      }
      coefficients( i, c ) = sum;
    }
  }
}

/// out = sum_{j >= first} coefficients(j, column) * src[j]
template< typename VECTOR >
void combine( std::vector< VECTOR * > const & src,
              DenseMatrix const & coefficients,
              integer const column,
              integer const first,
              VECTOR & out )
{
  out.zero();
  for( integer j = first; j < LvArray::integerConversion< integer >( src.size() ); ++j )
  {
    out.axpy( coefficients( j, column ), *src[j] );
  }
}

} // namespace

template< typename VECTOR >
EigenSolverResult LobpcgEigenSolver< VECTOR >::solve( Problem const & problem,
                                                      Vector const & prototype,
                                                      std::vector< Vector > & modes ) const
{
  GEOS_ERROR_IF( problem.preconditioner == nullptr && problem.shiftedInverse == nullptr,
                 "The LOBPCG eigensolver requires a preconditioner (an approximation of (K - shift M)^{-1})" );
  LinearOperator< VECTOR > const & preconditioner = problem.preconditioner != nullptr ? *problem.preconditioner
                                                                                      : *problem.shiftedInverse;

  EigenSolverParameters const & params = this->m_params;
  GEOS_ERROR_IF_LT_MSG( params.numEigenvalues, 1, "The number of requested eigenvalues must be positive" );

  Stopwatch watch;

  integer const nev = params.numEigenvalues;
  integer const n = std::max( nev, params.subspaceSize );
  globalIndex const dimension = problem.numUnknowns > 0 ? problem.numUnknowns : prototype.globalSize();
  GEOS_ERROR_IF( static_cast< globalIndex >( n ) >= dimension,
                 GEOS_FMT( "The problem has {} unknowns, which is too small for {} eigenpairs with LOBPCG", dimension, n ) );

  auto makeBlock = [&]()
  {
    std::vector< Vector > block;
    block.reserve( static_cast< size_t >( n ) );
    for( integer i = 0; i < n; ++i )
    {
      block.push_back( Base::makeVector( prototype ) );
    }
    return block;
  };

  // Iterates, search directions and preconditioned residuals, with their images by K and M
  std::vector< Vector > X = makeBlock();
  std::vector< Vector > KX = makeBlock();
  std::vector< Vector > MX = makeBlock();
  std::vector< Vector > Xn = makeBlock();
  std::vector< Vector > KXn = makeBlock();
  std::vector< Vector > MXn = makeBlock();
  std::vector< Vector > P = makeBlock();
  std::vector< Vector > KP = makeBlock();
  std::vector< Vector > MP = makeBlock();
  std::vector< Vector > Pn = makeBlock();
  std::vector< Vector > KPn = makeBlock();
  std::vector< Vector > MPn = makeBlock();
  std::vector< Vector > W = makeBlock();
  std::vector< Vector > KW = makeBlock();
  std::vector< Vector > MW = makeBlock();

  integer numOperatorApplications = 0;

  // Starting vectors: random, passed through the preconditioner so that they live in its range
  for( integer i = 0; i < n; ++i )
  {
    X[i].rand( static_cast< unsigned >( params.seed ) + 7919u * static_cast< unsigned >( i ) );
    problem.mass.apply( X[i], MX[i] );
    X[i].zero();
    preconditioner.apply( MX[i], X[i] );
    ++numOperatorApplications;
  }

  DenseMatrix coefficients;
  std::vector< real64 > theta;
  {
    std::vector< Vector * > S;
    std::vector< Vector * > KS;
    std::vector< Vector * > MS;
    for( integer i = 0; i < n; ++i )
    {
      problem.mass.apply( X[i], MX[i] );
      problem.stiffness.apply( X[i], KX[i] );
      S.push_back( &X[i] );
      KS.push_back( &KX[i] );
      MS.push_back( &MX[i] );
    }
    rayleighRitz( S, KS, MS, n, coefficients, theta );
    for( integer c = 0; c < n; ++c )
    {
      combine( S, coefficients, c, 0, Xn[c] );
      problem.mass.apply( Xn[c], MXn[c] );
      problem.stiffness.apply( Xn[c], KXn[c] );
    }
    std::swap( X, Xn );
    std::swap( KX, KXn );
    std::swap( MX, MXn );
  }

  integer iteration = 0;
  integer numConverged = 0;
  bool havePrevious = false;
  std::vector< real64 > errorEstimate( static_cast< size_t >( n ), 0.0 );

  for(;; )
  {
    // Preconditioned residuals W = T ( K x - theta M x ). Their M-norm estimates ||(K - sigma M)^{-1} r||_M, the
    // error measure of the shift-and-invert Krylov solver, and is not limited by the rounding noise of K x
    // that dominates the plain residual of the near-zero eigenvalues.
    std::vector< integer > active;
    numConverged = 0;
    for( integer i = 0; i < n; ++i )
    {
      KW[i].copy( KX[i] );
      KW[i].axpy( -theta[i], MX[i] );
      W[i].zero();
      preconditioner.apply( KW[i], W[i] );
      ++numOperatorApplications;
      problem.mass.apply( W[i], MW[i] );
      real64 const xNorm = std::sqrt( std::max( X[i].dot( MX[i] ), 0.0 ) );
      real64 const wNorm = std::sqrt( std::max( W[i].dot( MW[i] ), 0.0 ) );
      errorEstimate[i] = xNorm > 0.0 ? wNorm / xNorm : wNorm;
      bool const converged = i < nev && errorEstimate[i] <= params.tolerance;
      if( converged )
      {
        ++numConverged;
      }
      else
      {
        active.push_back( i );
      }
    }

    GEOS_LOG_RANK_0_IF( params.logLevel >= 1,
                        GEOS_FMT( "  LOBPCG: iteration {:4}, converged {:3}/{}, max error estimate {:.2e}, operator applications {}",
                                  iteration, numConverged, nev,
                                  *std::max_element( errorEstimate.begin(), errorEstimate.begin() + nev ),
                                  numOperatorApplications ) );

    if( numConverged == nev || iteration >= params.maxIterations )
    {
      break;
    }
    ++iteration;

    // Search space: [X, preconditioned residuals, previous directions] of the active columns
    std::vector< Vector * > S;
    std::vector< Vector * > KS;
    std::vector< Vector * > MS;
    for( integer i = 0; i < n; ++i )
    {
      S.push_back( &X[i] );
      KS.push_back( &KX[i] );
      MS.push_back( &MX[i] );
    }
    for( integer const i : active )
    {
      problem.stiffness.apply( W[i], KW[i] );
      S.push_back( &W[i] );
      KS.push_back( &KW[i] );
      MS.push_back( &MW[i] );
    }
    if( havePrevious )
    {
      for( integer const i : active )
      {
        S.push_back( &P[i] );
        KS.push_back( &KP[i] );
        MS.push_back( &MP[i] );
      }
    }

    rayleighRitz( S, KS, MS, n, coefficients, theta );

    for( integer c = 0; c < n; ++c )
    {
      // The images of the iterates are recomputed instead of being updated with the coefficients: the
      // recurrences drift when the search directions become nearly dependent, which limits the accuracy.
      combine( S, coefficients, c, 0, Xn[c] );
      problem.mass.apply( Xn[c], MXn[c] );
      problem.stiffness.apply( Xn[c], KXn[c] );
      combine( S, coefficients, c, n, Pn[c] );
      combine( KS, coefficients, c, n, KPn[c] );
      combine( MS, coefficients, c, n, MPn[c] );
    }
    std::swap( X, Xn );
    std::swap( KX, KXn );
    std::swap( MX, MXn );
    std::swap( P, Pn );
    std::swap( KP, KPn );
    std::swap( MP, MPn );
    havePrevious = true;
  }

  modes.clear();
  modes.reserve( static_cast< size_t >( nev ) );
  for( integer i = 0; i < nev; ++i )
  {
    modes.push_back( std::move( X[i] ) );
  }

  EigenSolverResult result;
  result.converged = ( numConverged == nev );
  result.numConverged = numConverged;
  result.numIterations = iteration;
  result.numOperatorApplications = numOperatorApplications;
  this->finalizeResult( problem, modes, result );
  result.solveTime = watch.elapsedTime();
  return result;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class LobpcgEigenSolver< TrilinosInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_HYPRE
template class LobpcgEigenSolver< HypreInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_PETSC
template class LobpcgEigenSolver< PetscInterface::ParallelVector >;
#endif

} // namespace geos

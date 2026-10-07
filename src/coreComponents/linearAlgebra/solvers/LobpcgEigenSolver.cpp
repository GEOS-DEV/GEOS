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
 * @return the dimension of the search space after removal of the dependent directions. If it is smaller than
 *         @p wanted, nothing is extracted and the caller needs to enlarge the basis.
 *
 * The basis is M-orthonormalized with SVQB: the Gram matrix is scaled to a unit diagonal, diagonalized, and the
 * directions with a negligible eigenvalue are dropped. The Ritz problem is then a standard symmetric eigenproblem.
 */
template< typename VECTOR >
integer rayleighRitz( stdVector< VECTOR * > const & S,
                      stdVector< VECTOR * > const & KS,
                      stdVector< VECTOR * > const & MS,
                      integer const wanted,
                      DenseMatrix & coefficients,
                      stdVector< real64 > & theta )
{
  integer const m = LvArray::integerConversion< integer >( S.size() );

  // Gram matrices of the search space. Each is computed by one batched device operation.
  array2d< real64 > products;
  DenseMatrix G( m, m );
  DenseMatrix A( m, m );
  stdVector< VECTOR const * > const basis = multiVectorOperations::constPointers( S );
  // Both Gram matrices are symmetric, so only their upper triangles are computed
  multiVectorOperations::dots( basis, multiVectorOperations::constPointers( MS ), products, true );
  for( integer j = 0; j < m; ++j )
  {
    for( integer i = 0; i < m; ++i )
    {
      G( i, j ) = 0.5 * ( products( i, j ) + products( j, i ) );
    }
  }
  multiVectorOperations::dots( basis, multiVectorOperations::constPointers( KS ), products, true );
  for( integer j = 0; j < m; ++j )
  {
    for( integer i = 0; i < m; ++i )
    {
      A( i, j ) = 0.5 * ( products( i, j ) + products( j, i ) );
    }
  }

  // Scale to a unit diagonal and diagonalize
  stdVector< real64 > d( m );
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

  stdVector< integer > kept;
  for( integer j = 0; j < m; ++j )
  {
    if( sigma[j] > svqbTolerance * sigma[m - 1] )
    {
      kept.push_back( j );
    }
  }
  integer const r = LvArray::integerConversion< integer >( kept.size() );
  if( r < wanted )
  {
    return r;
  }

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
  return r;
}

/// out = sum_{j >= first} coefficients(j, column) * src[j], in one fused device kernel
template< typename VECTOR >
void combine( stdVector< VECTOR * > const & src,
              DenseMatrix const & coefficients,
              integer const column,
              integer const first,
              VECTOR & out )
{
  stdVector< VECTOR const * > vectors;
  stdVector< real64 > weights;
  for( integer j = first; j < LvArray::integerConversion< integer >( src.size() ); ++j )
  {
    vectors.push_back( src[j] );
    weights.push_back( coefficients( j, column ) );
  }
  multiVectorOperations::combine( vectors, weights, out, false );
}

} // namespace

template< typename VECTOR >
EigenSolverResult LobpcgEigenSolver< VECTOR >::solve( Problem const & problem,
                                                      Vector const & prototype,
                                                      stdVector< Vector > & modes ) const
{
  GEOS_ERROR_IF( problem.preconditioner == nullptr && problem.shiftedInverse == nullptr,
                 "The LOBPCG eigensolver requires a preconditioner (an approximation of (K - shift M)^{-1})" );
  LinearOperator< VECTOR > const & preconditioner = problem.preconditioner != nullptr ? *problem.preconditioner
                                                                                      : *problem.shiftedInverse;

  EigenSolverParameters const & params = this->m_params;
  GEOS_ERROR_IF_LT_MSG( params.numEigenvalues, 1, "The number of requested eigenvalues must be positive" );

  Stopwatch watch;

  // Deflate the known eigenvectors: the iteration runs in their M-orthogonal complement
  ConstraintSpace constraints = this->makeConstraintSpace( problem, prototype );
  integer const numConstraints = constraints.size();
  integer const nev = params.numEigenvalues - numConstraints;
  if( nev <= 0 )
  {
    EigenSolverResult trivial = this->returnConstraintsOnly( problem, constraints, modes );
    trivial.solveTime = watch.elapsedTime();
    return trivial;
  }
  integer const n = std::max( nev, params.subspaceSize - numConstraints );
  globalIndex const dimension = ( problem.numUnknowns > 0 ? problem.numUnknowns : prototype.globalSize() ) - numConstraints;
  GEOS_ERROR_IF( static_cast< globalIndex >( n ) >= dimension,
                 GEOS_FMT( "The problem has {} unknowns, which is too small for {} eigenpairs with LOBPCG", dimension, n ) );

  auto makeBlock = [&]()
  {
    stdVector< Vector > block;
    block.reserve( static_cast< size_t >( n ) );
    for( integer i = 0; i < n; ++i )
    {
      block.push_back( Base::makeVector( prototype ) );
    }
    return block;
  };

  // Iterates, search directions and preconditioned residuals, with their images by K and M
  stdVector< Vector > X = makeBlock();
  stdVector< Vector > KX = makeBlock();
  stdVector< Vector > MX = makeBlock();
  stdVector< Vector > Xn = makeBlock();
  stdVector< Vector > KXn = makeBlock();
  stdVector< Vector > MXn = makeBlock();
  stdVector< Vector > P = makeBlock();
  stdVector< Vector > KP = makeBlock();
  stdVector< Vector > MP = makeBlock();
  stdVector< Vector > Pn = makeBlock();
  stdVector< Vector > KPn = makeBlock();
  stdVector< Vector > MPn = makeBlock();
  stdVector< Vector > W = makeBlock();
  stdVector< Vector > KW = makeBlock();
  stdVector< Vector > MW = makeBlock();

  DenseMatrix coefficients;
  stdVector< real64 > theta;
  integer numOperatorApplications = 0;
  unsigned randomCount = 0;

  // Random vector that vanishes on the constrained unknowns
  auto randomize = [&]( Vector & v )
  {
    v.rand( static_cast< unsigned >( params.seed ) + 7919u * randomCount++ );
    if( problem.freeMask != nullptr )
    {
      v.pointwiseProduct( *problem.freeMask );
    }
  };

  // Rayleigh-Ritz extraction. If the search space has fewer than n independent directions, it is enlarged with
  // random vectors, which are kept in padding for the duration of the call.
  stdVector< std::unique_ptr< Vector > > padding;
  auto extract = [&]( stdVector< Vector * > & S, stdVector< Vector * > & KS, stdVector< Vector * > & MS )
  {
    padding.clear();
    for( int attempt = 0; attempt < 10; ++attempt )
    {
      integer const rank = rayleighRitz( S, KS, MS, n, coefficients, theta );
      if( rank >= n )
      {
        return;
      }
      for( integer t = rank; t <= n; ++t )
      {
        padding.push_back( std::make_unique< Vector >( Base::makeVector( prototype ) ) );
        padding.push_back( std::make_unique< Vector >( Base::makeVector( prototype ) ) );
        padding.push_back( std::make_unique< Vector >( Base::makeVector( prototype ) ) );
        Vector & z = *padding[padding.size() - 3];
        randomize( z );
        constraints.project( z );
        problem.stiffness.apply( z, *padding[padding.size() - 2] );
        problem.mass.apply( z, *padding[padding.size() - 1] );
        S.push_back( &z );
        KS.push_back( padding[padding.size() - 2].get() );
        MS.push_back( padding[padding.size() - 1].get() );
      }
    }
    GEOS_ERROR( "LOBPCG: the search space has too few independent directions" );
  };

  // Starting vectors: random, with a unit-diagonal-scaled orthonormalization performed by the first extraction
  for( integer i = 0; i < n; ++i )
  {
    randomize( X[i] );
    constraints.project( X[i] );
  }

  {
    stdVector< Vector * > S;
    stdVector< Vector * > KS;
    stdVector< Vector * > MS;
    for( integer i = 0; i < n; ++i )
    {
      problem.mass.apply( X[i], MX[i] );
      problem.stiffness.apply( X[i], KX[i] );
      S.push_back( &X[i] );
      KS.push_back( &KX[i] );
      MS.push_back( &MX[i] );
    }
    extract( S, KS, MS );
    for( integer c = 0; c < n; ++c )
    {
      combine( S, coefficients, c, 0, Xn[c] );
      constraints.project( Xn[c] );
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
  stdVector< real64 > errorEstimate( static_cast< size_t >( n ), 0.0 );

  // A column that has converged is locked: it stays in the Rayleigh-Ritz space and does not make search
  // directions. Testing a column costs one preconditioner application, so a locked column is only tested every
  // lockedTestInterval iterations, which saves most of the cycles when the columns converge at different
  // iterations. The Rayleigh-Ritz steps mix the columns, so a locked column can drift: a locked column that fails
  // its test is unlocked. All the locked columns are also tested before the solver stops.
  integer constexpr lockedTestInterval = 5;
  stdVector< integer > locked( static_cast< size_t >( n ), false );

  // Preconditioned residual W = T ( K x - theta M x ) of column i. Its M-norm estimates ||(K - sigma M)^{-1} r||_M,
  // the error measure of the shift-and-invert Krylov solver, and is not limited by the rounding noise of K x
  // that dominates the plain residual of the near-zero eigenvalues.
  auto const testColumn = [&]( integer const i )
  {
    KW[i].copy( KX[i] );
    KW[i].axpy( -theta[i], MX[i] );
    W[i].zero();
    preconditioner.apply( KW[i], W[i] );
    ++numOperatorApplications;
    constraints.project( W[i] );
    problem.mass.apply( W[i], MW[i] );
    real64 const xNorm = std::sqrt( std::max( X[i].dot( MX[i] ), 0.0 ) );
    real64 const wNorm = std::sqrt( std::max( W[i].dot( MW[i] ), 0.0 ) );
    errorEstimate[i] = xNorm > 0.0 ? wNorm / xNorm : wNorm;
    return i < nev && errorEstimate[i] <= params.tolerance;
  };

  for(;; )
  {
    stdVector< integer > active;
    stdVector< integer > tested( static_cast< size_t >( n ), false );
    numConverged = 0;
    for( integer i = 0; i < n; ++i )
    {
      if( locked[i] && iteration % lockedTestInterval != 0 )
      {
        ++numConverged;
        continue;
      }
      tested[i] = true;
      if( testColumn( i ) )
      {
        locked[i] = true;
        ++numConverged;
      }
      else
      {
        locked[i] = false;
        active.push_back( i );
      }
    }

    // Before stopping, test the locked columns that were not tested in this iteration
    if( numConverged == nev || iteration >= params.maxIterations )
    {
      bool unlocked = false;
      for( integer i = 0; i < nev; ++i )
      {
        if( locked[i] && !tested[i] )
        {
          if( !testColumn( i ) )
          {
            locked[i] = false;
            --numConverged;
            active.push_back( i );
            unlocked = true;
          }
        }
      }
      if( unlocked )
      {
        std::sort( active.begin(), active.end() );
      }
    }

    GEOS_LOG_RANK_0_IF( params.logLevel >= 1,
                        GEOS_FMT( "  LOBPCG: iteration {:4}, converged {:3}/{}, max error estimate {:.2e}, operator applications {}",
                                  iteration, numConverged, nev,
                                  *std::max_element( errorEstimate.begin(), errorEstimate.begin() + nev ),
                                  numOperatorApplications ) );

    if( params.logLevel >= 2 )
    {
      string line;
      for( integer i = 0; i < nev; ++i )
      {
        line += GEOS_FMT( " {:.1e}", errorEstimate[i] );
      }
      GEOS_LOG_RANK_0( GEOS_FMT( "    error estimates:{}", line ) );
    }

    if( numConverged == nev || iteration >= params.maxIterations )
    {
      break;
    }
    ++iteration;

    // Search space: [X, preconditioned residuals, previous directions] of the active columns
    stdVector< Vector * > S;
    stdVector< Vector * > KS;
    stdVector< Vector * > MS;
    for( integer i = 0; i < n; ++i )
    {
      S.push_back( &X[i] );
      KS.push_back( &KX[i] );
      MS.push_back( &MX[i] );
    }

    // Separate small search directions from the O(1) iterate block before SVQB.
    // Otherwise near-parallel P and X columns amplify roundoff in the projected
    // pencil and impose a residual floor on consistent-mass elasticity problems.
    stdVector< Vector * > directions, massDirections;
    for( integer const i : active )
    {
      directions.push_back( &W[i] );
      massDirections.push_back( &MW[i] );
      if( havePrevious )
      {
        directions.push_back( &P[i] );
        massDirections.push_back( &MP[i] );
      }
    }
    stdVector< Vector const * > const iterates = multiVectorOperations::constPointers( S );
    stdVector< Vector const * > const massIterates = multiVectorOperations::constPointers( MS );
    // The images by M of the directions are known when entering (W: computed by the test of the column, P: refreshed
    // after the last extraction). Each pass updates them with the same combination as the directions, and one
    // application at the end removes the drift.
    for( integer pass = 0; pass < 2; ++pass )
    {
      array2d< real64 > products;
      multiVectorOperations::dots( iterates, multiVectorOperations::constPointers( massDirections ), products );
      for( size_t j = 0; j < directions.size(); ++j )
      {
        stdVector< real64 > weights( static_cast< size_t >( n ) );
        for( integer i = 0; i < n; ++i )
          weights[i] = -products( i, LvArray::integerConversion< localIndex >( j ) );
        multiVectorOperations::combine( iterates, weights, *directions[j], true );
        multiVectorOperations::combine( massIterates, weights, *massDirections[j], true );
      }
    }
    for( size_t j = 0; j < directions.size(); ++j )
      problem.mass.apply( *directions[j], *massDirections[j] );

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
        problem.stiffness.apply( P[i], KP[i] );
        S.push_back( &P[i] );
        KS.push_back( &KP[i] );
        MS.push_back( &MP[i] );
      }
    }

    extract( S, KS, MS );

    for( integer c = 0; c < n; ++c )
    {
      // The images of the iterates are recomputed instead of being updated with the coefficients: the
      // recurrences drift when the search directions become nearly dependent, which limits the accuracy.
      combine( S, coefficients, c, 0, Xn[c] );
      constraints.project( Xn[c] );
      problem.mass.apply( Xn[c], MXn[c] );
      problem.stiffness.apply( Xn[c], KXn[c] );
      combine( S, coefficients, c, n, Pn[c] );
      constraints.project( Pn[c] );
      // Search directions can be formed by cancellation of almost parallel vectors.
      // Updating their images by the same recurrence then ceases to represent K P and M P,
      // corrupting the next projected pencil. Refresh them just like the iterate images.
      problem.stiffness.apply( Pn[c], KPn[c] );
      problem.mass.apply( Pn[c], MPn[c] );
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
  modes.reserve( static_cast< size_t >( nev + numConstraints ) );
  for( integer i = 0; i < numConstraints; ++i )
  {
    modes.push_back( std::move( constraints.vectors[i] ) );
  }
  for( integer i = 0; i < nev; ++i )
  {
    modes.push_back( std::move( X[i] ) );
  }

  EigenSolverResult result;
  result.converged = ( numConverged == nev );
  result.numConverged = numConverged + numConstraints;
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

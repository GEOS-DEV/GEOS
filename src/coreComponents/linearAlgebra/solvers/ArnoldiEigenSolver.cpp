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
 * @file ArnoldiEigenSolver.cpp
 */

#include "ArnoldiEigenSolver.hpp"

#include "common/Stopwatch.hpp"
#include "denseLinearAlgebra/interfaces/blaslapack/BlasLapackLA.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <numeric>

namespace geos
{

namespace
{

/// Ritz pairs of the projected operator, ordered by decreasing modulus of the Ritz value
struct RitzPairs
{
  /// Dimension of the Rayleigh-Ritz problem
  integer m = 0;
  /// Ritz values theta of the operator T = (K - sigma M)^{-1} M
  std::vector< real64 > theta;
  /// Coefficients of the Ritz vectors in the Krylov basis (column i is the vector of theta[i])
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > Y;
  /// Residual estimates ||T x_i - theta_i x_i||_M
  std::vector< real64 > rho;
  /// Coefficients of the residuals in the (not yet expanded) next block of the basis: coupling[c + b*i]
  std::vector< real64 > coupling;
};

/**
 * @brief State of a block Krylov-Schur iteration.
 *
 * The basis V[0..m+b-1] is M-orthonormal; V[0..m-1] spans the Krylov space, and V[m..m+b-1] is the next block,
 * not yet expanded. The (m+b) x m matrix H holds the coefficients of T V[j] on the basis, so that
 * T V_m = V_{m+b} H. The Rayleigh-Ritz matrix is the symmetrized leading m x m block of H.
 */
template< typename VECTOR >
class KrylovSchurState
{
public:

  using Problem = GeneralizedEigenProblem< VECTOR >;

  KrylovSchurState( EigenSolverParameters const & params,
                    Problem const & problem,
                    VECTOR const & prototype,
                    integer const ncv,
                    integer const blockSize ):
    m_params( params ),
    m_problem( problem ),
    m_ncv( ncv ),
    m_b( blockSize ),
    m_ld( ncv + blockSize ),
    m_H( static_cast< size_t >( m_ld ) * static_cast< size_t >( ncv ), 0.0 )
  {
    size_t const numBasis = static_cast< size_t >( ncv + blockSize );
    m_V.reserve( numBasis );
    m_MV.reserve( numBasis );
    for( size_t i = 0; i < numBasis; ++i )
    {
      m_V.push_back( makeVector( prototype ) );
      m_MV.push_back( makeVector( prototype ) );
    }
    m_X.reserve( static_cast< size_t >( ncv ) );
    for( integer i = 0; i < ncv; ++i )
    {
      m_X.push_back( makeVector( prototype ) );
    }
  }

  /// @return the number of applications of the shift-and-invert operator so far
  integer numOperatorApplications() const { return m_numOperatorApplications; }

  /// Fill the first block of the basis with random vectors in the range of T
  void initialize()
  {
    for( integer c = 0; c < m_b; ++c )
    {
      randomVector( c );
    }
  }

  /// Append one block of b vectors to the basis. @p m is the current number of columns of H.
  void extendBlock( integer const m )
  {
    std::vector< real64 > coefficients( static_cast< size_t >( m_ld ) );
    for( integer c = 0; c < m_b; ++c )
    {
      integer const column = m + c;
      integer const index = m + m_b + c;

      applyOperator( m_MV[column], m_V[index] );

      std::fill( coefficients.begin(), coefficients.end(), 0.0 );
      real64 const norm = orthonormalize( index, coefficients.data() );
      for( integer i = 0; i < index; ++i )
      {
        h( i, column ) = coefficients[i];
      }
      if( norm > 0.0 )
      {
        h( index, column ) = norm;
      }
      else
      {
        // The Krylov space is invariant (or numerically so): continue with a fresh direction
        h( index, column ) = 0.0;
        randomVector( index );
      }
    }
  }

  /// Compute the Ritz pairs of the Rayleigh-Ritz problem of dimension @p m
  RitzPairs ritz( integer const m ) const
  {
    array2d< real64, MatrixLayout::COL_MAJOR_PERM > S( m, m );
    array2d< real64, MatrixLayout::COL_MAJOR_PERM > Yraw( m, m );
    array1d< real64 > lambda( m );
    for( integer j = 0; j < m; ++j )
    {
      for( integer i = 0; i <= j; ++i )
      {
        S( i, j ) = h( i, j );
        S( j, i ) = h( i, j );
      }
    }
    BlasLapackLA::matrixSymmetricEigen( S.toSliceConst(), lambda.toSlice(), Yraw.toSlice() );

    std::vector< integer > order( static_cast< size_t >( m ) );
    std::iota( order.begin(), order.end(), 0 );
    std::stable_sort( order.begin(), order.end(), [&]( integer const a, integer const b )
    {
      return std::fabs( lambda[a] ) > std::fabs( lambda[b] );
    } );

    RitzPairs r;
    r.m = m;
    r.theta.resize( m );
    r.rho.resize( m );
    r.coupling.assign( static_cast< size_t >( m_b ) * static_cast< size_t >( m ), 0.0 );
    r.Y.resize( m, m );
    for( integer i = 0; i < m; ++i )
    {
      integer const o = order[i];
      r.theta[i] = lambda[o];
      for( integer j = 0; j < m; ++j )
      {
        r.Y( j, i ) = Yraw( j, o );
      }
      real64 norm2 = 0.0;
      for( integer c = 0; c < m_b; ++c )
      {
        real64 s = 0.0;
        for( integer j = 0; j < m; ++j )
        {
          s += h( m + c, j ) * r.Y( j, i );
        }
        r.coupling[c + m_b * i] = s;
        norm2 += s * s;
      }
      r.rho[i] = std::sqrt( norm2 );
    }
    return r;
  }

  /**
   * @brief Number of leading Ritz pairs satisfying rho <= tol * |theta|.
   * @param r the Ritz pairs
   * @param wanted number of leading pairs to examine
   * @return the number of converged pairs
   */
  integer countConverged( RitzPairs const & r, integer const wanted ) const
  {
    real64 const floor = std::pow( DBL_EPSILON, 2.0 / 3.0 );
    integer n = 0;
    for( integer i = 0; i < wanted; ++i )
    {
      if( r.rho[i] <= m_params.tolerance * std::max( std::fabs( r.theta[i] ), floor ) )
      {
        ++n;
      }
    }
    return n;
  }

  /**
   * @brief Krylov-Schur restart.
   * @param k number of Ritz vectors kept
   * @param r Ritz pairs of the current basis of dimension r.m
   * @param injectRandomBlock if true, the next block is random (orthogonal to the kept vectors) instead of the
   *        residual direction of the current basis
   *
   * Afterwards the basis has k + b vectors and H is arrow-shaped.
   */
  void restart( integer const k, RitzPairs const & r, bool const injectRandomBlock )
  {
    integer const m = r.m;
    ritzVectors( k, r );

    if( !injectRandomBlock )
    {
      for( integer c = 0; c < m_b; ++c )
      {
        m_V[k + c].copy( m_V[m + c] );
        m_MV[k + c].copy( m_MV[m + c] );
      }
    }
    for( integer i = 0; i < k; ++i )
    {
      m_V[i].copy( m_X[i] );
      m_problem.mass.apply( m_V[i], m_MV[i] );
    }

    std::fill( m_H.begin(), m_H.end(), 0.0 );
    for( integer i = 0; i < k; ++i )
    {
      h( i, i ) = r.theta[i];
    }
    if( injectRandomBlock )
    {
      for( integer c = 0; c < m_b; ++c )
      {
        randomVector( k + c );
      }
    }
    else
    {
      for( integer c = 0; c < m_b; ++c )
      {
        for( integer i = 0; i < k; ++i )
        {
          h( k + c, i ) = r.coupling[c + m_b * i];
          h( i, k + c ) = r.coupling[c + m_b * i];
        }
      }
    }
  }

  /// Move the first @p count Ritz vectors to @p modes
  void extractModes( integer const count, RitzPairs const & r, std::vector< VECTOR > & modes )
  {
    ritzVectors( count, r );
    modes.clear();
    modes.reserve( static_cast< size_t >( count ) );
    for( integer i = 0; i < count; ++i )
    {
      modes.push_back( std::move( m_X[i] ) );
    }
  }

private:

  static VECTOR makeVector( VECTOR const & prototype )
  {
    VECTOR v;
    v.create( prototype.localSize(), prototype.comm() );
    v.zero();
    return v;
  }

  real64 & h( integer const i, integer const j ) { return m_H[ i + static_cast< size_t >( j ) * m_ld ]; }
  real64 h( integer const i, integer const j ) const { return m_H[ i + static_cast< size_t >( j ) * m_ld ]; }

  /// dst = (K - sigma M)^{-1} src
  void applyOperator( VECTOR const & src, VECTOR & dst )
  {
    dst.zero();
    m_problem.shiftedInverse->apply( src, dst );
    ++m_numOperatorApplications;
  }

  /// m_X[0..count-1] = V[0..m-1] Y[:, 0..count-1]
  void ritzVectors( integer const count, RitzPairs const & r )
  {
    for( integer i = 0; i < count; ++i )
    {
      m_X[i].zero();
      for( integer j = 0; j < r.m; ++j )
      {
        m_X[i].axpy( r.Y( j, i ), m_V[j] );
      }
    }
  }

  /**
   * @brief M-orthonormalize V[index] against V[0..index-1] (modified Gram-Schmidt, two passes).
   * @param index index of the vector to process
   * @param coefficients if non-null, receives the Gram-Schmidt coefficients (accumulated over the passes)
   * @return the M-norm of the vector after orthogonalization, or zero if it is numerically in the span of
   *         the previous vectors (the vector is then left unnormalized)
   */
  real64 orthonormalize( integer const index, real64 * const coefficients )
  {
    VECTOR & w = m_V[index];
    VECTOR & mw = m_MV[index];

    m_problem.mass.apply( w, mw );
    real64 const norm0 = std::sqrt( std::max( w.dot( mw ), 0.0 ) );
    if( norm0 <= 0.0 )
    {
      return 0.0;
    }

    for( int pass = 0; pass < 2; ++pass )
    {
      for( integer i = 0; i < index; ++i )
      {
        real64 const c = m_MV[i].dot( w );
        w.axpy( -c, m_V[i] );
        if( coefficients != nullptr )
        {
          coefficients[i] += c;
        }
      }
    }

    m_problem.mass.apply( w, mw );
    real64 const norm = std::sqrt( std::max( w.dot( mw ), 0.0 ) );
    if( norm <= breakdownTolerance * norm0 )
    {
      return 0.0;
    }
    w.scale( 1.0 / norm );
    mw.scale( 1.0 / norm );
    return norm;
  }

  /// Fill V[index] with a random vector in the range of T, M-orthonormal to V[0..index-1]
  void randomVector( integer const index )
  {
    for( int attempt = 0; attempt < 5; ++attempt )
    {
      m_V[index].rand( static_cast< unsigned >( m_params.seed ) + 7919u * static_cast< unsigned >( m_randomCount++ ) );
      m_problem.mass.apply( m_V[index], m_MV[index] );
      applyOperator( m_MV[index], m_V[index] );
      if( orthonormalize( index, nullptr ) > 0.0 )
      {
        return;
      }
    }
    GEOS_ERROR( "Eigensolver could not generate a new starting vector: the problem is too small for the "
                "requested subspace size, or M is singular on the whole space." );
  }

  /// Relative drop of the norm below which a vector is considered to be in the span of the basis
  static constexpr real64 breakdownTolerance = 1.0e-8;

  EigenSolverParameters const & m_params;
  Problem const & m_problem;
  integer const m_ncv;
  integer const m_b;
  integer const m_ld;
  std::vector< real64 > m_H;
  std::vector< VECTOR > m_V;
  std::vector< VECTOR > m_MV;
  std::vector< VECTOR > m_X;
  integer m_numOperatorApplications = 0;
  integer m_randomCount = 0;
};

} // namespace

template< typename VECTOR >
EigenSolverResult ArnoldiEigenSolver< VECTOR >::solve( Problem const & problem,
                                                       Vector const & prototype,
                                                       std::vector< Vector > & modes ) const
{
  GEOS_ERROR_IF( problem.shiftedInverse == nullptr,
                 "The Arnoldi eigensolver requires the shift-and-invert operator (K - shift M)^{-1}" );

  EigenSolverParameters const & params = this->m_params;
  GEOS_ERROR_IF_LT_MSG( params.numEigenvalues, 1, "The number of requested eigenvalues must be positive" );
  GEOS_ERROR_IF_LT_MSG( params.blockSize, 1, "The eigensolver block size must be positive" );

  Stopwatch watch;

  integer const nev = params.numEigenvalues;
  integer const b = params.blockSize;

  // Basis dimension: a multiple of the block size, large enough to hold the wanted vectors and a block, and
  // small enough for the basis to be linearly independent
  integer ncv = params.subspaceSize > 0 ? params.subspaceSize : std::max( 2 * nev, nev + 2 * b );
  ncv = std::max( ncv, nev + b );
  ncv = ( ncv + b - 1 ) / b * b;
  globalIndex const maxBasis = problem.numUnknowns > 0 ? problem.numUnknowns : prototype.globalSize();
  if( static_cast< globalIndex >( ncv + b ) > maxBasis )
  {
    ncv = LvArray::integerConversion< integer >( ( maxBasis - b ) / b * b );
  }
  GEOS_ERROR_IF( ncv < nev + b,
                 GEOS_FMT( "The problem has {} unknowns, which is too small for {} eigenpairs with block size {}",
                           maxBasis, nev, b ) );

  KrylovSchurState< Vector > state( params, problem, prototype, ncv, b );
  state.initialize();

  integer constexpr maxCompletenessChecks = 3;
  integer m = 0;
  integer restarts = 0;
  integer numChecks = 0;
  integer numConverged = 0;
  bool verifying = false;
  RitzPairs ritz;

  for(;; )
  {
    bool converged = false;
    bool haveRitz = false;
    while( m + b <= ncv )
    {
      state.extendBlock( m );
      m += b;
      haveRitz = false;
      if( !verifying && m >= nev )
      {
        ritz = state.ritz( m );
        haveRitz = true;
        if( state.countConverged( ritz, nev ) == nev )
        {
          converged = true;
          break;
        }
      }
    }
    if( !haveRitz )
    {
      ritz = state.ritz( m );
    }
    numConverged = state.countConverged( ritz, nev );
    converged = ( numConverged == nev );

    GEOS_LOG_RANK_0_IF( params.logLevel >= 1,
                        GEOS_FMT( "  Arnoldi {}: restart {:3}, basis {:3}, converged {:3}/{}, operator applications {}",
                                  verifying ? "check" : "     ", restarts, m, numConverged, nev,
                                  state.numOperatorApplications() ) );

    if( converged )
    {
      if( verifying || params.completenessCheck == 0 || numChecks >= maxCompletenessChecks )
      {
        break;
      }
      // Look for missed copies of repeated eigenvalues: keep the converged vectors and expand a random block
      // that is orthogonal to them
      ++numChecks;
      verifying = true;
      state.restart( nev, ritz, true );
      m = nev;
      continue;
    }

    if( restarts >= params.maxIterations )
    {
      break;
    }
    ++restarts;
    verifying = false;
    integer const keep = nev + std::min( numConverged, ( ncv - b - nev ) / 2 );
    state.restart( keep, ritz, false );
    m = keep;
  }

  state.extractModes( nev, ritz, modes );

  EigenSolverResult result;
  result.converged = ( numConverged == nev );
  result.numConverged = numConverged;
  result.numIterations = restarts;
  result.numOperatorApplications = state.numOperatorApplications();
  this->finalizeResult( problem, modes, result );
  result.solveTime = watch.elapsedTime();
  return result;
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOS_USE_TRILINOS
template class ArnoldiEigenSolver< TrilinosInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_HYPRE
template class ArnoldiEigenSolver< HypreInterface::ParallelVector >;
#endif

#ifdef GEOS_USE_PETSC
template class ArnoldiEigenSolver< PetscInterface::ParallelVector >;
#endif

} // namespace geos

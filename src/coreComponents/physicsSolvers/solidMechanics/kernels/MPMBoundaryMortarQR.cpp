/* SPDX-License-Identifier: LGPL-2.1-only
 *
 * Numerical backends for the existing boundary-admissible minimum-mass
 * projection. No GEOS or MPI types: compile this translation unit directly in
 * the standalone tests. In GEOS, CMake supplies the optional library definitions.
 */
#include "MPMBoundaryMortarQR.hpp"

#ifdef GEOS_MPM_MORTAR_LAPACK
// Match GEOS' selected Fortran naming convention. This backend intentionally
// supports LP64 (32-bit LAPACK integers), not ILP64 libraries with the same names.
#ifdef GEOS_MPM_MORTAR_FORTRAN_NO_UNDERSCORE
#define MPM_LAPACK(name) name
#else
#define MPM_LAPACK(name) name ## _
#endif
extern "C"
{
void MPM_LAPACK(dgeqp3)( int const *, int const *, double *, int const *, int *,
                         double *, double *, int const *, int * );
void MPM_LAPACK(dormqr)( char const *, char const *, int const *, int const *, int const *,
                         double const *, int const *, double const *, double *,
                         int const *, double *, int const *, int * );
void MPM_LAPACK(dtrtrs)( char const *, char const *, char const *, int const *, int const *,
                         double const *, int const *, double *, int const *, int * );
}
#endif

namespace geos
{
namespace mpm
{
namespace boundaryMortar
{
namespace
{

void validateDimensions( std::vector< double > const & bt, Index n, Index m,
                         std::vector< double > const & b, double tolerance )
{
  require( m == 0 || n <= std::numeric_limits< Index >::max() / m,
           "Boundary mortar matrix dimensions overflow." );
  require( b.size() == m && bt.size() == n*m, "Invalid rectangular QR dimensions." );
  require( std::isfinite( tolerance ) && tolerance > 0.0 && tolerance < 1.0,
           "Invalid boundary mortar QR rank threshold." );
  for( double value : bt ) require( std::isfinite( value ), "Nonfinite QR coefficient." );
  for( double value : b ) require( std::isfinite( value ), "Nonfinite QR right-hand side." );
}

// project() normalizes each column of B^T. Sum-of-squares is safe for those
// inputs and permits vectorization. Rank decisions use refreshed, not stale norms.
double trailingNorm( std::vector< double > const & a, Index n, Index column, Index start )
{
  double sum = 0.0;
  for( Index i = start; i < n; ++i ) sum += a[column*n+i]*a[column*n+i];
  return std::sqrt( sum );
}

#ifdef GEOS_MPM_MORTAR_LAPACK
int workspaceSize( double query, Index minimum )
{
  require( std::isfinite( query ) && query >= 1.0 &&
           query <= double( std::numeric_limits< int >::max() ) &&
           minimum <= Index( std::numeric_limits< int >::max() ),
           "LAPACK workspace exceeds the supported LP64 limit." );
  return std::max( int( std::ceil( query ) ), int( minimum ) );
}
#endif

} // namespace

bool lapackQRAvailable()
{
#ifdef GEOS_MPM_MORTAR_LAPACK
  return true;
#else
  return false;
#endif
}

QRBackend selectQRBackend( std::string const & name )
{
  if( name == "auto" ) return lapackQRAvailable() ? QRBackend::Lapack : QRBackend::Native;
  if( name == "reference" ) return QRBackend::Reference;
  if( name == "native" ) return QRBackend::Native;
  if( name == "lapack" )
  {
    require( lapackQRAvailable(),
             "weakInterfaceTraceQRBackend=lapack requires GEOS_ENABLE_MPM_MORTAR_LAPACK=ON." );
    return QRBackend::Lapack;
  }
  throw std::runtime_error( "weakInterfaceTraceQRBackend must be auto, reference, native, or lapack." );
}

char const * qrBackendName( QRBackend backend )
{
  switch( backend )
  {
    case QRBackend::Reference: return "reference";
    case QRBackend::Native: return "native";
    case QRBackend::Lapack: return "lapack";
  }
  throw std::runtime_error( "Invalid boundary mortar QR backend." );
}

MinimumNormSolver minimumNormSolver( QRBackend backend )
{
  switch( backend )
  {
    case QRBackend::Reference: return minimumNorm;
    case QRBackend::Native: return minimumNormNative;
    case QRBackend::Lapack:
      require( lapackQRAvailable(), "LAPACK mortar backend is unavailable in this build." );
      return minimumNormLapack;
  }
  throw std::runtime_error( "Invalid boundary mortar QR backend." );
}

QRResult minimumNormNative( std::vector< double > & bt, Index const n, Index const m,
                            std::vector< double > const & b, double const rankTolerance )
{
  validateDimensions( bt, n, m, b, rankTolerance );
  QRResult out;
  out.correction.assign( n, 0.0 ); out.dual.assign( m, 0.0 );
  Index const limit = std::min( n, m );
  std::vector< Index > permutation( m );
  std::iota( permutation.begin(), permutation.end(), Index( 0 ) );
  std::vector< double > tau( limit ), norms( m ), exactNorms( m );
  for( Index j = 0; j < m; ++j ) norms[j] = exactNorms[j] = trailingNorm( bt, n, j, 0 );
  double firstPivot = 0.0;
  double const refreshTolerance = std::sqrt( std::numeric_limits< double >::epsilon() );

  for( Index k = 0; k < limit; ++k )
  {
    Index pivot = k;
    for( Index j = k+1; j < m; ++j ) if( norms[j] > norms[pivot] ) pivot = j;

    // Recheck every remaining norm before dropping a direction. Partial-norm
    // estimates are an acceleration, never authority to relax the rank test.
    if( k > 0 && !( norms[pivot] > rankTolerance*firstPivot ) )
    {
      for( Index j = k; j < m; ++j ) norms[j] = exactNorms[j] = trailingNorm( bt, n, j, k );
      pivot = k;
      for( Index j = k+1; j < m; ++j ) if( norms[j] > norms[pivot] ) pivot = j;
    }
    double const pivotNorm = trailingNorm( bt, n, pivot, k );
    if( k == 0 ) firstPivot = pivotNorm;
    if( !( pivotNorm > rankTolerance*firstPivot ) || !( pivotNorm > 0.0 ) )
    {
      // A near-threshold estimate may have selected an underestimated column.
      // Refresh and reconsider, rather than declaring numerical rank early.
      for( Index j = k; j < m; ++j ) norms[j] = exactNorms[j] = trailingNorm( bt, n, j, k );
      pivot = k;
      for( Index j = k+1; j < m; ++j ) if( norms[j] > norms[pivot] ) pivot = j;
      if( !( norms[pivot] > rankTolerance*firstPivot ) || !( norms[pivot] > 0.0 ) ) break;
    }
    double const bestNorm = trailingNorm( bt, n, pivot, k );
    if( pivot != k )
    {
      for( Index i = 0; i < n; ++i ) std::swap( bt[k*n+i], bt[pivot*n+i] );
      std::swap( permutation[k], permutation[pivot] );
      std::swap( norms[k], norms[pivot] );
      std::swap( exactNorms[k], exactNorms[pivot] );
    }
    double const x0 = bt[k*n+k];
    double const alpha = -std::copysign( bestNorm, x0 );
    double const denominator = x0-alpha;
    tau[k] = (alpha-x0)/alpha;
    for( Index i = k+1; i < n; ++i ) bt[k*n+i] /= denominator;
    bt[k*n+k] = alpha;
    for( Index j = k+1; j < m; ++j )
    {
      double dot = bt[j*n+k];
      for( Index i = k+1; i < n; ++i ) dot += bt[k*n+i]*bt[j*n+i];
      dot *= tau[k];
      bt[j*n+k] -= dot;
      for( Index i = k+1; i < n; ++i ) bt[j*n+i] -= bt[k*n+i]*dot;

      // Safeguarded norm downdating (LAPACK Working Note 176 / DLAQP2).
      // Refresh when cancellation makes the estimate unreliable.
      if( norms[j] > 0.0 )
      {
        double const ratio = std::abs( bt[j*n+k] )/norms[j];
        double const remaining = std::max( 0.0, (1.0+ratio)*(1.0-ratio) );
        double const history = norms[j]/exactNorms[j];
        if( remaining*history*history <= refreshTolerance )
          norms[j] = exactNorms[j] = trailingNorm( bt, n, j, k+1 );
        else
          norms[j] *= std::sqrt( remaining );
      }
    }
    ++out.rank;
    out.smallestRelativePivot = bestNorm/firstPivot;
  }

  // Same minimum-norm and dual reconstruction as the reference implementation.
  // Pivot order may differ at floating-point ties; the physical correction and
  // the original-row residual, not bitwise multipliers, are the contract.
  Index const r = out.rank;
  std::vector< double > z( r, 0.0 );
  for( Index i = 0; i < r; ++i )
  {
    double rhs = b[permutation[i]];
    for( Index j = 0; j < i; ++j ) rhs -= bt[i*n+j]*z[j];
    z[i] = rhs/bt[i*n+i]; out.correction[i] = z[i];
  }
  for( Index kk = r; kk > 0; --kk )
  {
    Index const k = kk-1;
    double dot = out.correction[k];
    for( Index i = k+1; i < n; ++i ) dot += bt[k*n+i]*out.correction[i];
    dot *= tau[k]; out.correction[k] -= dot;
    for( Index i = k+1; i < n; ++i ) out.correction[i] -= bt[k*n+i]*dot;
  }
  std::vector< double > eta( r, 0.0 );
  for( Index ii = r; ii > 0; --ii )
  {
    Index const i = ii-1;
    double rhs = z[i];
    for( Index j = i+1; j < r; ++j ) rhs -= bt[j*n+i]*eta[j];
    eta[i] = rhs/bt[i*n+i]; out.dual[permutation[i]] = eta[i];
  }
  return out;
}

QRResult minimumNormLapack( std::vector< double > & bt, Index const n, Index const m,
                            std::vector< double > const & b, double const rankTolerance )
{
#ifndef GEOS_MPM_MORTAR_LAPACK
  (void)bt; (void)n; (void)m; (void)b; (void)rankTolerance;
  throw std::runtime_error( "LAPACK mortar backend was not compiled." );
#else
  static_assert( sizeof( int ) == 4, "The optional mortar LAPACK backend requires LP64 integers." );
  validateDimensions( bt, n, m, b, rankTolerance );
  require( n <= Index( std::numeric_limits< int >::max() ) &&
           m <= Index( (std::numeric_limits< int >::max()-1)/3 ),
           "Boundary mortar QR exceeds LAPACK LP64 dimension/workspace limits." );
  QRResult out;
  out.correction.assign( n, 0.0 ); out.dual.assign( m, 0.0 );
  if( n == 0 || m == 0 ) return out;
  int const nn = int( n ), mm = int( m ), lda = nn, one = 1;
  std::vector< int > permutation( m, 0 ); // all columns free to pivot
  std::vector< double > tau( std::min( n, m ) );
  int info = 0, lwork = -1;
  double query = 0.0;
  MPM_LAPACK(dgeqp3)( &nn, &mm, bt.data(), &lda, permutation.data(), tau.data(), &query, &lwork, &info );
  require( info == 0, "DGEQP3 workspace query failed." );
  lwork = workspaceSize( query, 3*m+1 );
  std::vector< double > work( Index( lwork ), 0.0 );
  MPM_LAPACK(dgeqp3)( &nn, &mm, bt.data(), &lda, permutation.data(), tau.data(), work.data(), &lwork, &info );
  require( info == 0, "DGEQP3 factorization failed." );
  double const firstPivot = std::abs( bt[0] );
  for( Index k = 0; k < tau.size(); ++k )
  {
    double const pivot = std::abs( bt[k*n+k] );
    require( std::isfinite( pivot ), "Nonfinite LAPACK QR pivot." );
    if( !( pivot > rankTolerance*firstPivot ) || !( pivot > 0.0 ) ) break;
    ++out.rank; out.smallestRelativePivot = pivot/firstPivot;
  }
  int const r = int( out.rank );
  if( r == 0 ) return out;
  std::vector< double > z( out.rank );
  for( int i = 0; i < r; ++i )
  {
    require( permutation[i] >= 1 && permutation[i] <= mm, "Invalid LAPACK QR permutation." );
    z[i] = b[permutation[i]-1];
  }
  char const upper = 'U', transpose = 'T', noTranspose = 'N', nonUnit = 'N', left = 'L';
  MPM_LAPACK(dtrtrs)( &upper, &transpose, &nonUnit, &r, &one, bt.data(), &lda, z.data(), &r, &info );
  require( info == 0, "QR transpose triangular solve failed." );
  std::copy( z.begin(), z.end(), out.correction.begin() );

  // Apply only the retained reflectors to [z;0]. Do not form Q explicitly.
  lwork = -1; query = 0.0;
  MPM_LAPACK(dormqr)( &left, &noTranspose, &nn, &one, &r, bt.data(), &lda, tau.data(),
                     out.correction.data(), &nn, &query, &lwork, &info );
  require( info == 0, "DORMQR workspace query failed." );
  lwork = workspaceSize( query, 1 );
  work.resize( Index( lwork ) );
  MPM_LAPACK(dormqr)( &left, &noTranspose, &nn, &one, &r, bt.data(), &lda, tau.data(),
                     out.correction.data(), &nn, work.data(), &lwork, &info );
  require( info == 0, "DORMQR application failed." );

  // Keep a valid dual for existing impulse/reaction accounting. Dependent-row
  // multipliers need not match another pivoting scheme; their primal impulse must.
  MPM_LAPACK(dtrtrs)( &upper, &noTranspose, &nonUnit, &r, &one, bt.data(), &lda, z.data(), &r, &info );
  require( info == 0, "QR dual triangular solve failed." );
  for( int i = 0; i < r; ++i ) out.dual[permutation[i]-1] = z[i];
  return out;
#endif
}

} // namespace boundaryMortar
} // namespace mpm
} // namespace geos

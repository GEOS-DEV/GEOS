/*
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Small CPU reference solver for a boundary-compatible MPM mortar projection.
 * No GEOS, MPI, or linear-algebra-library types are used here so the algebra can
 * be tested independently. This is intentionally a bounded, dense reference
 * solve, NOT a distributed/scalable replacement for GEOS' linear solver stack.
 */
#ifndef GEOS_MPM_BOUNDARY_MORTAR_HPP
#define GEOS_MPM_BOUNDARY_MORTAR_HPP

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace geos
{
namespace mpm
{
namespace boundaryMortar
{

using Index = std::size_t;
constexpr Index invalidIndex = std::numeric_limits< Index >::max();

struct Options
{
  double relativeTolerance = 1.0e-10;
  double absoluteTolerance = 1.0e-12; // velocity units, for the final average jump
  double rankTolerance = 1.0e-12;     // relative QR pivot threshold, NOT compliance
  double minimumRetainedFraction = 1.0e-12;
  Index maximumMatrixEntries = 16000000; // 128 MB for the factor itself
  bool supportFallback = true;
};

// A constraint is the integral of a velocity jump. Each coefficient multiplies
// (v[a] - v[b]); both sides retain exactly the same quadrature coefficient.
struct PairCoefficient
{
  Index a = 0;
  Index b = 0;
  double value = 0.0;
};
struct Row
{
  std::vector< PairCoefficient > coefficients;
  double measure = 0.0;
};

// Either v[a] = sign*v[b] + offset, or v[a] = offset (b == invalidIndex).
// Normal symmetry/moving-buffer reflection has sign=-1; nodal fallback sign=+1.
struct Relation
{
  Index a = 0;
  Index b = invalidIndex;
  double sign = 1.0;
  double offset = 0.0;
};

inline void require( bool const condition, char const * message )
{
  if( !condition ) throw std::runtime_error( message );
}

/** Signed affine disjoint sets: v_i = sign_i*q_root + offset_i.
 *
 * This is the sparse Z map for the six grid-aligned GEOS faces. One raw degree
 * of freedom has at most one free root, so Z^T M Z is a sum of masses, not a
 * dense matrix. Fixed values and closed reflection cycles are handled exactly.
 * A contradictory prescription is an error, never an averaged compromise.
 */
class AffineMap
{
public:
  explicit AffineMap( Index const n, double const tolerance )
    : m_parent( n ), m_size( n, 1 ), m_sign( n, 1.0 ), m_offset( n, 0.0 ),
      m_fixed( n, false ), m_value( n, 0.0 ), m_tolerance( tolerance )
  {
    std::iota( m_parent.begin(), m_parent.end(), Index( 0 ) );
  }

  Index root( Index const i )
  {
    require( i < m_parent.size(), "Boundary dependency index is out of range." );
    if( m_parent[i] != i )
    {
      Index const p = m_parent[i];
      m_parent[i] = root( p );
      m_offset[i] += m_sign[i] * m_offset[p];
      m_sign[i] *= m_sign[p];
    }
    return m_parent[i];
  }

  void fix( Index const i, double const value )
  {
    require( std::isfinite( value ), "Nonfinite prescribed boundary velocity." );
    Index const r = root( i );
    setRootValue( r, m_sign[i] * ( value - m_offset[i] ) );
  }

  void relate( Index const a, Index const b, double const sign, double const offset )
  {
    bool const positiveSign = close( sign, 1.0 );
    bool const negativeSign = close( sign, -1.0 );
    require( ( positiveSign || negativeSign ) && std::isfinite( offset ),
             "Boundary relation must have a finite offset and a unit sign." );
    double const relationSign = positiveSign ? 1.0 : -1.0;
    Index const ra = root( a );
    Index const rb = root( b );
    double const s = m_sign[a] * relationSign * m_sign[b];
    double const d = m_sign[a] * ( relationSign * m_offset[b] + offset - m_offset[a] );
    // q_ra = s*q_rb + d.
    if( ra == rb )
    {
      if( s > 0.0 )
        require( close( d, 0.0 ), "Incompatible boundary/fallback dependency cycle." );
      else
        setRootValue( ra, 0.5 * d );
      return;
    }
    if( m_size[ra] <= m_size[rb] )
      attach( ra, rb, s, d );
    else
      attach( rb, ra, s, -s * d );
  }

  void apply( Relation const & relation )
  {
    if( relation.b == invalidIndex ) fix( relation.a, relation.offset );
    else relate( relation.a, relation.b, relation.sign, relation.offset );
  }

  double sign( Index const i ) { root( i ); return m_sign[i]; }
  double offset( Index const i ) { root( i ); return m_offset[i]; }
  bool fixed( Index const i ) { return m_fixed[root( i )]; }
  double fixedValue( Index const i ) { return m_value[root( i )]; }

private:
  bool close( double const a, double const b ) const
  {
    return std::abs( a - b ) <= m_tolerance * std::max( { 1.0, std::abs( a ), std::abs( b ) } );
  }
  void setRootValue( Index const r, double const value )
  {
    if( m_fixed[r] )
      require( close( m_value[r], value ), "Incompatible prescribed boundary velocities." );
    m_fixed[r] = true;
    m_value[r] = value;
  }
  void attach( Index const child, Index const parent, double const sign, double const offset )
  {
    if( m_fixed[child] ) setRootValue( parent, sign * ( m_value[child] - offset ) );
    m_parent[child] = parent;
    m_sign[child] = sign;
    m_offset[child] = offset;
    m_size[parent] += m_size[child];
  }
  std::vector< Index > m_parent, m_size;
  std::vector< double > m_sign, m_offset;
  std::vector< bool > m_fixed;
  std::vector< double > m_value;
  double m_tolerance;
};

struct ReducedState
{
  std::vector< Index > root, column;
  std::vector< double > sign, offset, mass, center, value;
  std::vector< bool > fixed;
  Index freeCount = 0;
};

inline ReducedState reduce( AffineMap & map, std::vector< double > const & mass,
                            std::vector< double > const & trial )
{
  Index const n = mass.size();
  ReducedState s;
  s.root.resize( n ); s.column.assign( n, invalidIndex );
  s.sign.resize( n ); s.offset.resize( n ); s.mass.assign( n, 0.0 );
  s.center.assign( n, 0.0 ); s.value.resize( n ); s.fixed.resize( n );
  for( Index i = 0; i < n; ++i )
  {
    require( mass[i] >= 0.0 && std::isfinite( mass[i] ) && std::isfinite( trial[i] ),
             "Invalid mass or trial velocity in the boundary mortar solve." );
    Index const r = map.root( i );
    s.root[i] = r; s.sign[i] = map.sign( i ); s.offset[i] = map.offset( i );
    s.fixed[i] = map.fixed( i );
    s.mass[r] += mass[i];
    s.center[r] += mass[i] * s.sign[i] * ( trial[i] - s.offset[i] );
  }
  for( Index r = 0; r < n; ++r )
  {
    if( s.root[r] != r ) continue;
    if( s.fixed[r] ) s.center[r] = map.fixedValue( r );
    else if( s.mass[r] > 0.0 )
    {
      s.center[r] /= s.mass[r];
      s.column[r] = s.freeCount++;
    }
    // Unsupported free roots never enter a trace row. Keep the representative
    // trial value rather than clearing an unrelated zero-mass velocity.
    else s.center[r] = trial[r];
  }
  for( Index i = 0; i < n; ++i )
    s.value[i] = s.sign[i] * s.center[s.root[i]] + s.offset[i];
  return s;
}

struct QRResult
{
  std::vector< double > correction, dual;
  Index rank = 0;
  double smallestRelativePivot = 0.0;
};

/** Solve min ||y|| subject to B*y=b with column-pivoted Householder QR of B^T.
 * B^T P = Q [R11 R12; 0 0]. For a compatible right-hand side the unique minimum
 * norm correction is y=Q1 R11^{-T} (P^T b)_1. The caller verifies ALL original
 * constraints, including dependent rows. A rank threshold never adds compliance.
 *
 * bt is column-major: bt[column_of_Bt*n + row_of_Bt]. It is overwritten by QR.
 */
inline QRResult minimumNorm( std::vector< double > & bt, Index const n, Index const m,
                             std::vector< double > const & b, double const rankTolerance )
{
  require( b.size() == m && bt.size() == n * m, "Invalid rectangular QR dimensions." );
  QRResult out;
  out.correction.assign( n, 0.0 ); out.dual.assign( m, 0.0 );
  std::vector< Index > permutation( m );
  std::iota( permutation.begin(), permutation.end(), Index( 0 ) );
  Index const limit = std::min( n, m );
  std::vector< double > tau( limit, 0.0 );
  double firstPivot = 0.0;
  for( Index k = 0; k < limit; ++k )
  {
    Index pivot = k;
    double bestNorm = -1.0;
    // Recompute norms rather than downdating: slower, but robust for this
    // verification-sized, strongly correlated mortar matrix.
    for( Index j = k; j < m; ++j )
    {
      // project() normalizes every column of B^T to unit norm before
      // factorization. Orthogonal transforms preserve that bound, so this
      // dot product cannot overflow; no expensive per-entry hypot is needed.
      double normSquared = 0.0;
      for( Index i = k; i < n; ++i ) normSquared += bt[j*n+i] * bt[j*n+i];
      double const norm = std::sqrt( normSquared );
      if( norm > bestNorm ) { bestNorm = norm; pivot = j; }
    }
    if( k == 0 ) firstPivot = bestNorm;
    if( !( bestNorm > rankTolerance * firstPivot ) || !( bestNorm > 0.0 ) ) break;
    if( pivot != k )
    {
      for( Index i = 0; i < n; ++i ) std::swap( bt[k*n+i], bt[pivot*n+i] );
      std::swap( permutation[k], permutation[pivot] );
    }
    double const x0 = bt[k*n+k];
    double const alpha = -std::copysign( bestNorm, x0 );
    double const denominator = x0 - alpha;
    tau[k] = ( alpha - x0 ) / alpha;
    for( Index i = k + 1; i < n; ++i ) bt[k*n+i] /= denominator;
    bt[k*n+k] = alpha;
    for( Index j = k + 1; j < m; ++j )
    {
      double dot = bt[j*n+k];
      for( Index i = k + 1; i < n; ++i ) dot += bt[k*n+i] * bt[j*n+i];
      dot *= tau[k];
      bt[j*n+k] -= dot;
      for( Index i = k + 1; i < n; ++i ) bt[j*n+i] -= bt[k*n+i] * dot;
    }
    ++out.rank;
    out.smallestRelativePivot = bestNorm / firstPivot;
  }
  Index const r = out.rank;
  std::vector< double > z( r, 0.0 );
  for( Index i = 0; i < r; ++i )
  {
    double rhs = b[permutation[i]];
    for( Index j = 0; j < i; ++j ) rhs -= bt[i*n+j] * z[j];
    z[i] = rhs / bt[i*n+i];
    out.correction[i] = z[i];
  }
  // Q = H0 H1 ...: apply Householder reflectors in reverse order to Q1*z.
  for( Index kk = r; kk > 0; --kk )
  {
    Index const k = kk - 1;
    double dot = out.correction[k];
    for( Index i = k + 1; i < n; ++i ) dot += bt[k*n+i] * out.correction[i];
    dot *= tau[k];
    out.correction[k] -= dot;
    for( Index i = k + 1; i < n; ++i ) out.correction[i] -= bt[k*n+i] * dot;
  }
  // One valid dual solution; multipliers need not be unique when rows repeat.
  std::vector< double > eta( r, 0.0 );
  for( Index ii = r; ii > 0; --ii )
  {
    Index const i = ii - 1;
    double rhs = z[i];
    for( Index j = i + 1; j < r; ++j ) rhs -= bt[j*n+i] * eta[j];
    eta[i] = rhs / bt[i*n+i];
    out.dual[permutation[i]] = eta[i];
  }
  return out;
}

struct Result
{
  std::vector< double > velocity, mortarImpulse, reactionImpulse;
  std::vector< double > preJump, postJump, retainedFraction;
  std::vector< int > accepted, fallbackDof;
  Index rank = 0, unknowns = 0, assembledRows = 0, rejectedRows = 0;
  double maximumJump = 0.0, relativeResidual = 0.0, boundaryResidual = 0.0;
  double smallestRelativePivot = 0.0;
};

/** Boundary reduction, common-support selection, and simultaneous mortar solve.
 *
 * mass contains owned physical nodal-field masses (duplicates have already been
 * removed). Zero-mass buffer entries may depend on mass-carrying donor entries.
 * Support is checked AFTER boundary reduction. New support fallback ties are
 * incorporated into that reduction; there is no post-solve velocity overwrite.
 */
inline Result project( std::vector< double > const & mass,
                       std::vector< double > const & trial,
                       std::vector< double > const & cutoff,
                       std::vector< Relation > const & boundary,
                       std::vector< std::pair< Index, Index > > const & fallback,
                       std::vector< Row > const & rows,
                       Options const & options )
{
  Index const n = mass.size();
  require( trial.size() == n && cutoff.size() == n, "Boundary mortar input size mismatch." );
  require( std::isfinite( options.relativeTolerance ) && std::isfinite( options.absoluteTolerance ) &&
           std::isfinite( options.rankTolerance ) &&
           options.relativeTolerance > 0.0 && options.absoluteTolerance > 0.0 &&
           options.rankTolerance > 0.0 && options.rankTolerance < 1.0 &&
           options.minimumRetainedFraction >= 0.0 && options.minimumRetainedFraction <= 1.0,
           "Invalid boundary mortar solver tolerances." );
  for( double const value : cutoff )
    require( value >= 0.0 && std::isfinite( value ), "Invalid boundary mortar mass cutoff." );
  AffineMap map( n, options.absoluteTolerance );
  for( Relation const & r : boundary ) map.apply( r );
  Result result;
  result.fallbackDof.assign( n, 0 );
  auto tie = [&]( Index const a, Index const b )
  {
    require( a < n && b < n, "Fallback index out of range." );
    map.relate( a, b, 1.0, 0.0 );
    result.fallbackDof[a] = result.fallbackDof[b] = 1;
  };
  for( auto const & pair : fallback )
  {
    require( pair.first < n && pair.second < n, "Fallback index out of range." );
    if( mass[pair.first] > cutoff[pair.first] && mass[pair.second] > cutoff[pair.second] )
      tie( pair.first, pair.second );
  }

  result.accepted.assign( rows.size(), 1 );
  result.retainedFraction.assign( rows.size(), 0.0 );
  std::vector< std::vector< PairCoefficient > > retained( rows.size() );
  ReducedState state;
  // At most one new row removal is required for progress. This also handles the
  // rare case where a fallback merges previously separate boundary donor sets.
  for( Index pass = 0; pass <= rows.size(); ++pass )
  {
    state = reduce( map, mass, trial );
    bool removed = false;
    for( Index a = 0; a < rows.size(); ++a )
    {
      if( !result.accepted[a] ) continue;
      Row const & row = rows[a];
      require( row.measure > 0.0 && std::isfinite( row.measure ), "Invalid mortar row measure." );
      double total = 0.0, kept = 0.0;
      retained[a].clear();
      for( PairCoefficient const & c : row.coefficients )
      {
        require( c.a < n && c.b < n && std::isfinite( c.value ),
                 "Invalid mortar coefficient." );
        // Preserve coefficients exactly. The GEOS partition-of-unity closure
        // can produce a roundoff-sized negative last B-spline coefficient.
        total += c.value;
        if( state.mass[state.root[c.a]] > cutoff[c.a] && state.mass[state.root[c.b]] > cutoff[c.b] )
        { retained[a].push_back( c ); kept += c.value; }
      }
      result.retainedFraction[a] = total > 0.0 ? kept / total : 0.0;
      if( kept <= 0.0 || result.retainedFraction[a] < options.minimumRetainedFraction )
      {
        result.accepted[a] = 0; ++result.rejectedRows; removed = true;
        if( options.supportFallback )
          for( PairCoefficient const & c : row.coefficients )
            if( mass[c.a] > cutoff[c.a] && mass[c.b] > cutoff[c.b] ) tie( c.a, c.b );
      }
    }
    if( !removed ) break;
    require( pass < rows.size(), "Support/fallback reduction failed to terminate." );
  }

  result.preJump.assign( rows.size(), 0.0 );
  result.postJump.assign( rows.size(), 0.0 );
  std::vector< std::map< Index, double > > reducedRows;
  std::vector< Index > sourceRow;
  std::vector< double > scale, rhs;
  double initialNorm = 0.0;
  for( Index a = 0; a < rows.size(); ++a )
  {
    if( !result.accepted[a] ) continue;
    std::map< Index, double > coefficients;
    double g = 0.0;
    for( PairCoefficient const & c : retained[a] )
    {
      g += c.value * ( state.value[c.a] - state.value[c.b] );
      for( int side = 0; side < 2; ++side )
      {
        Index const i = side == 0 ? c.a : c.b;
        Index const r = state.root[i];
        if( state.column[r] != invalidIndex )
          coefficients[state.column[r]] += ( side == 0 ? 1.0 : -1.0 ) * c.value *
                                           state.sign[i] / std::sqrt( state.mass[r] );
      }
    }
    result.preJump[a] = g / rows[a].measure;
    initialNorm = std::hypot( initialNorm, result.preJump[a] );
    double norm = 0.0;
    for( auto const & c : coefficients ) norm = std::hypot( norm, c.second );
    if( !( norm > 0.0 ) )
    {
      require( std::abs( result.preJump[a] ) <= options.absoluteTolerance,
               "Mortar row conflicts with fully prescribed boundary/fallback velocities." );
      continue; // redundant but satisfied; checked again below
    }
    sourceRow.push_back( a ); scale.push_back( 1.0 / norm ); rhs.push_back( -g / norm );
    reducedRows.emplace_back( std::move( coefficients ) );
  }
  // Remove free roots untouched by any accepted mortar row before allocating QR.
  // Their mass-minimizing boundary projection is still retained in state.value.
  std::set< Index > usedColumns;
  for( auto const & row : reducedRows )
    for( auto const & c : row ) if( std::abs( c.second ) > 0.0 ) usedColumns.insert( c.first );
  std::map< Index, Index > compact;
  for( Index c : usedColumns ) compact.emplace( c, compact.size() );
  Index const nc = compact.size(), nr = reducedRows.size();
  require( nr == 0 || nc <= options.maximumMatrixEntries / nr,
           "Boundary mortar QR exceeds weakInterfaceTraceQRMaxEntries; reduce the reference problem or use a scalable solver." );
  std::vector< double > bt( nc * nr, 0.0 );
  for( Index a = 0; a < nr; ++a )
    for( auto const & c : reducedRows[a] )
      if( std::abs( c.second ) > 0.0 ) bt[a*nc+compact.at( c.first )] = scale[a] * c.second;
  QRResult qr = minimumNorm( bt, nc, nr, rhs, options.rankTolerance );
  result.rank = qr.rank; result.unknowns = nc; result.assembledRows = nr;
  result.smallestRelativePivot = qr.smallestRelativePivot;
  for( Index r = 0; r < n; ++r )
  {
    if( state.root[r] != r || state.column[r] == invalidIndex ) continue;
    auto const found = compact.find( state.column[r] );
    if( found != compact.end() ) state.center[r] += qr.correction[found->second] / std::sqrt( state.mass[r] );
  }
  result.velocity.resize( n ); result.mortarImpulse.assign( n, 0.0 );
  result.reactionImpulse.resize( n );
  for( Index i = 0; i < n; ++i )
    result.velocity[i] = state.sign[i] * state.center[state.root[i]] + state.offset[i];
  for( Index a = 0; a < nr; ++a )
  {
    double const impulse = qr.dual[a] * scale[a];
    for( PairCoefficient const & c : retained[sourceRow[a]] )
    { result.mortarImpulse[c.a] += c.value * impulse; result.mortarImpulse[c.b] -= c.value * impulse; }
  }
  double finalNorm = 0.0;
  for( Index a = 0; a < rows.size(); ++a )
  {
    if( !result.accepted[a] ) continue;
    double g = 0.0;
    for( PairCoefficient const & c : retained[a] )
      g += c.value * ( result.velocity[c.a] - result.velocity[c.b] );
    result.postJump[a] = g / rows[a].measure;
    finalNorm = std::hypot( finalNorm, result.postJump[a] );
    result.maximumJump = std::max( result.maximumJump, std::abs( result.postJump[a] ) );
  }
  result.relativeResidual = initialNorm > 0.0 ? finalNorm / initialNorm : finalNorm;
  double const rms = finalNorm / std::sqrt( double( std::max( Index( 1 ), rows.size() - result.rejectedRows ) ) );
  double const initialRms = initialNorm / std::sqrt( double( std::max( Index( 1 ), rows.size() - result.rejectedRows ) ) );
  require( rms <= options.absoluteTolerance + options.relativeTolerance * initialRms &&
           result.maximumJump <= options.absoluteTolerance + options.relativeTolerance *
             ( result.preJump.empty() ? 0.0 : std::abs( *std::max_element( result.preJump.begin(), result.preJump.end(),
                 []( double const a, double const b ){ return std::abs( a ) < std::abs( b ); } ) ) ),
           "Boundary mortar QR did not satisfy the original accepted constraints; inspect rank, support, and tolerances." );
  for( Relation const & r : boundary )
  {
    double const residual = result.velocity[r.a] - r.offset -
      ( r.b == invalidIndex ? 0.0 : r.sign * result.velocity[r.b] );
    result.boundaryResidual = std::max( result.boundaryResidual, std::abs( residual ) );
  }
  require( result.boundaryResidual <= 8.0 * options.absoluteTolerance,
           "Boundary dependency reconstruction lost the prescribed kinematics." );
  for( Index i = 0; i < n; ++i )
  {
    require( std::isfinite( result.velocity[i] ), "Nonfinite boundary mortar solution." );
    // This includes both external BC reactions and internal fallback reactions.
    // The latter cancel when summed over the fields tied at a node.
    result.reactionImpulse[i] = mass[i] * ( result.velocity[i] - trial[i] ) - result.mortarImpulse[i];
  }
  return result;
}

/** Global-index description of the existing, grid-aligned GEOS BC surfaces.
 * The physical faces are one cell inside the nonperiodic padded grid.
 * Interface orientation does not enter these boundary relations.
 */
using LatticeKey = std::array< long long, 3 >;
struct Faces
{
  std::array< long long, 3 > cells{{ 0, 0, 0 }};
  std::array< double, 3 > low{{ 0, 0, 0 }}, high{{ 0, 0, 0 }};
  std::array< bool, 3 > periodic{{ false, false, false }};
  std::array< int, 6 > kind{{ 0, 0, 0, 0, 0, 0 }}; // 0 outflow, 1 symmetry, 2 Moving
  std::array< std::array< bool, 3 >, 6 > controlled{};
  std::array< std::array< double, 3 >, 6 > target{};
  int dimensions = 3;

  bool inside( std::array< double, 3 > const & x ) const
  {
    for( int d = 0; d < dimensions; ++d )
      if( !periodic[d] && ( x[d] < low[d] || x[d] > high[d] ) ) return false;
    return true;
  }
  int boundaryMask( std::array< double, 3 > const & x, double const tolerance ) const
  {
    int mask = 0;
    for( int d = 0; d < dimensions; ++d ) if( !periodic[d] )
    {
      if( std::abs( x[d] - low[d] ) <= tolerance ) mask |= 1 << ( 2*d );
      if( std::abs( x[d] - high[d] ) <= tolerance ) mask |= 1 << ( 2*d+1 );
    }
    return mask;
  }
};

// A donor alone is not sufficient: include every buffer entry tied to it so all
// physical mapped mass in Z^T M Z is present, even when that entry has c_aI=0.
// This closure is global-index based and therefore independent of MPI slicing.
inline void closeBoundaryDependencies( std::set< LatticeKey > & keys, Faces const & faces )
{
  bool changed = true;
  while( changed )
  {
    changed = false;
    std::vector< LatticeKey > additions;
    for( LatticeKey const & key : keys )
      for( int face = 0; face < 2 * faces.dimensions; ++face )
      {
        int const d = face / 2;
        if( faces.periodic[d] || faces.kind[face] == 0 ) continue;
        long long const ghost = face % 2 == 0 ? 0 : faces.cells[d];
        long long const donor = face % 2 == 0 ? 2 : faces.cells[d] - 2;
        if( key[d] == ghost || key[d] == donor )
        {
          LatticeKey reflected = key;
          reflected[d] = key[d] == ghost ? donor : ghost;
          additions.push_back( reflected );
        }
      }
    for( LatticeKey const & key : additions ) changed = keys.insert( key ).second || changed;
  }
}

inline std::vector< Relation > makeBoundaryRelations( std::vector< LatticeKey > const & keys,
                                                       Index const fields, int const component,
                                                       Faces const & faces )
{
  std::map< LatticeKey, Index > lookup;
  for( Index i = 0; i < keys.size(); ++i ) lookup.emplace( keys[i], i );
  std::vector< Relation > relations;
  for( Index i = 0; i < keys.size(); ++i )
  {
    std::array< int, 6 > priority{};
    int highest = 0;
    for( int face = 0; face < 2 * faces.dimensions; ++face )
    {
      int const d = face / 2;
      if( faces.periodic[d] || !faces.controlled[face][component] ) continue;
      long long const boundary = face % 2 == 0 ? 1 : faces.cells[d] - 1;
      long long const buffer = face % 2 == 0 ? 0 : faces.cells[d];
      // Same precedence as GEOS' staged normal/transverse boundary targets.
      if( keys[i][d] == boundary ) priority[face] = component == d ? 4 : 2;
      if( keys[i][d] == buffer ) priority[face] = component == d ? 3 : 1;
      highest = std::max( highest, priority[face] );
    }
    for( int face = 0; face < 2 * faces.dimensions; ++face )
    {
      if( priority[face] == 0 || priority[face] != highest ) continue;
      int const d = face / 2;
      bool const fixed = priority[face] == 4 || priority[face] == 2;
      Index donorIndex = invalidIndex;
      if( !fixed )
      {
        LatticeKey donor = keys[i];
        donor[d] = face % 2 == 0 ? 2 : faces.cells[d] - 2;
        auto const found = lookup.find( donor );
        require( found != lookup.end(), "Missing MPI donor in boundary dependency closure." );
        donorIndex = found->second;
      }
      for( Index field = 0; field < fields; ++field )
        relations.push_back( Relation{ i*fields+field,
          fixed ? invalidIndex : donorIndex*fields+field,
          fixed ? 1.0 : -1.0, ( fixed ? 1.0 : 2.0 ) * faces.target[face][component] } );
    }
  }
  return relations;
}

struct SurfaceSample
{
  std::array< double, 3 > point{}, normal{};
};

// Only symmetry is reflected geometrically. A Moving face does not imply a
// mirror-symmetric material geometry. Images are fit samples, never particles,
// volume/mass contributions, or integration facets outside the physical domain.
inline std::vector< SurfaceSample > symmetryImages( SurfaceSample const & sample,
                                                    Faces const & faces,
                                                    double const radius, double const tolerance )
{
  std::vector< SurfaceSample > images{ sample };
  for( int d = 0; d < faces.dimensions; ++d )
  {
    Index const count = images.size();
    for( int side = 0; side < 2; ++side )
    {
      int const face = 2*d + side;
      if( faces.kind[face] != 1 || faces.periodic[d] ) continue;
      double const plane = side == 0 ? faces.low[d] : faces.high[d];
      for( Index j = 0; j < count; ++j )
      {
        double const distance = images[j].point[d] - plane;
        if( std::abs( distance ) <= tolerance || std::abs( distance ) > radius ) continue;
        SurfaceSample mirrored = images[j];
        mirrored.point[d] = 2.0 * plane - mirrored.point[d];
        mirrored.normal[d] = -mirrored.normal[d];
        images.push_back( mirrored );
      }
    }
  }
  return images;
}

} // namespace boundaryMortar
} // namespace mpm
} // namespace geos
#endif

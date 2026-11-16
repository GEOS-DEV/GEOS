/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_MPM_THERMAL_NUMERICS_HPP_
#define GEOS_MPM_THERMAL_NUMERICS_HPP_

// Storage-independent host reference operators, shared by the solver and tests.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace geos
{
namespace mpmThermal
{
using Vector = std::array< double, 3 >;
inline double dot( Vector const & a, Vector const & b )
{ return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
inline Vector cross( Vector const & a, Vector const & b )
{ return {{ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] }}; }
inline bool normalize( Vector & n )
{
  double const length = std::sqrt( dot( n, n ) );
  if( !( length > 0.0 ) || !std::isfinite( length ) ) return false;
  for( double & x : n ) x /= length;
  return true;
}

struct FourierAssembly
{
  std::vector< double > power;
  std::vector< double > absoluteRowBound;
};

/** Particle quadrature of the weak Fourier operator and a Gershgorin bound. */
inline FourierAssembly assembleFourierParticle( double const volume, double const conductivity,
                                                std::vector< Vector > const & gradients,
                                                std::vector< double > const & temperature )
{
  if( gradients.size() != temperature.size() || !( volume > 0.0 ) ||
      !std::isfinite( volume ) || conductivity < 0.0 || !std::isfinite( conductivity ) )
    throw std::invalid_argument( "Invalid Fourier particle quadrature" );
  Vector gradT = {{ 0.0, 0.0, 0.0 }};
  double normSum = 0.0;
  for( std::size_t i = 0; i < gradients.size(); ++i )
  {
    for( int k = 0; k < 3; ++k ) gradT[k] += temperature[i]*gradients[i][k];
    normSum += std::sqrt( dot( gradients[i], gradients[i] ) );
  }
  FourierAssembly result;
  for( auto const & gradient : gradients )
  {
    result.power.push_back( -volume*conductivity*dot( gradient, gradT ) );
    result.absoluteRowBound.push_back( volume*conductivity*std::sqrt( dot( gradient, gradient ) )*normSum );
  }
  return result;
}

inline double gatherSpecificEnergy( double const specificHeat, double const weight,
                                     double const nodalEnergy, double const nodalCapacity )
{ return specificHeat*weight*nodalEnergy/nodalCapacity; }

class UnionFind
{
public:
  explicit UnionFind( std::size_t const n ) : m_parent( n )
  { std::iota( m_parent.begin(), m_parent.end(), std::size_t( 0 ) ); }
  std::size_t root( std::size_t i )
  {
    while( i != m_parent[i] ) { m_parent[i] = m_parent[m_parent[i]]; i = m_parent[i]; }
    return i;
  }
  void join( std::size_t a, std::size_t b )
  {
    a = root( a ); b = root( b );
    if( a != b ) m_parent[std::max( a, b )] = std::min( a, b );
  }
private:
  std::vector< std::size_t > m_parent;
};

struct ThermalNode
{
  double capacity = 0.0;
  double temperature = 0.0;
  double bulkEnergy = 0.0;
  bool prescribed = false;
  double prescribedTemperature = 0.0;
};

struct ThermalEdge
{
  std::size_t a;
  std::size_t b;
  double conductance;
  bool perfect = false;
};

/** Small dense SPD solve. Inputs are scaled to avoid an unnecessary 1/dt. */
inline std::vector< double > solveSPD( std::vector< double > a, std::vector< double > rhs )
{
  std::size_t const n = rhs.size();
  for( std::size_t i = 0; i < n; ++i )
  {
    for( std::size_t j = 0; j <= i; ++j )
    {
      double v = a[i*n+j];
      for( std::size_t k = 0; k < j; ++k ) v -= a[i*n+k]*a[j*n+k];
      if( i == j )
      {
        if( !( v > 0.0 ) || !std::isfinite( v ) )
          throw std::runtime_error( "Thermal contact system is not numerically positive definite; reduce the thermal substep" );
        a[i*n+j] = std::sqrt( v );
      }
      else a[i*n+j] = v/a[j*n+j];
    }
  }
  for( std::size_t i = 0; i < n; ++i )
  {
    for( std::size_t j = 0; j < i; ++j ) rhs[i] -= a[i*n+j]*rhs[j];
    rhs[i] /= a[i*n+i];
  }
  for( std::size_t ii = n; ii > 0; --ii )
  {
    std::size_t const i = ii-1;
    for( std::size_t j = i+1; j < n; ++j ) rhs[i] -= a[j*n+i]*rhs[j];
    rhs[i] /= a[i*n+i];
  }
  return rhs;
}

/** Backward Euler exchange on the complete eligible contact graph.
 * Return physical energy increments, including heat entering prescribed nodes.
 * The particle thermostat supplies the equal/opposite reservoir reaction later.
 * Perfect edges are eliminated exactly; a zero resistance is never divided by.
 */
inline std::vector< double > exchangeGraph( double const dt,
                                           std::vector< ThermalNode > const & nodes,
                                           std::vector< ThermalEdge > const & edges )
{
  if( !( dt > 0.0 ) || !std::isfinite( dt ) ) throw std::invalid_argument( "Invalid thermal substep" );
  std::size_t const n = nodes.size();
  UnionFind sets( n );
  for( auto const & edge : edges )
  {
    if( edge.a >= n || edge.b >= n || edge.a == edge.b ||
        !std::isfinite( edge.conductance ) || edge.conductance < 0.0 )
      throw std::invalid_argument( "Invalid thermal contact edge" );
    if( edge.perfect ) sets.join( edge.a, edge.b );
  }
  std::vector< std::size_t > roots, index( n );
  for( std::size_t i = 0; i < n; ++i )
  {
    if( !( nodes[i].capacity > 0.0 ) || !std::isfinite( nodes[i].capacity ) ||
        !std::isfinite( nodes[i].temperature ) || !std::isfinite( nodes[i].bulkEnergy ) )
      throw std::invalid_argument( "Invalid thermal nodal state" );
    auto const root = sets.root( i );
    auto it = std::find( roots.begin(), roots.end(), root );
    if( it == roots.end() ) { roots.push_back( root ); index[i] = roots.size()-1; }
    else index[i] = static_cast< std::size_t >( it-roots.begin() );
  }
  std::size_t const m = roots.size();
  std::vector< double > capacity( m, 0.0 ), rhs( m, 0.0 ), fixedTemperature( m, 0.0 );
  std::vector< bool > fixed( m, false );
  for( std::size_t i = 0; i < n; ++i )
  {
    auto const r = index[i];
    capacity[r] += nodes[i].capacity;
    rhs[r] += nodes[i].capacity*nodes[i].temperature + nodes[i].bulkEnergy;
    if( nodes[i].prescribed )
    {
      double const t = nodes[i].prescribedTemperature;
      if( !std::isfinite( t ) || t < 0.0 || ( fixed[r] &&
          std::abs( fixedTemperature[r]-t ) > 1.0e-10*std::max( 1.0, std::abs( t ) ) ) )
        throw std::invalid_argument( "Conflicting temperatures on a perfectly tied thermal interface" );
      fixed[r] = true; fixedTemperature[r] = t;
    }
  }
  std::vector< double > matrix( m*m, 0.0 );
  for( std::size_t i = 0; i < m; ++i ) matrix[i*m+i] = capacity[i];
  for( auto const & edge : edges )
  {
    auto const a = index[edge.a], b = index[edge.b];
    if( a == b || edge.perfect ) continue;
    double const g = dt*edge.conductance;
    if( !std::isfinite( g ) ) throw std::overflow_error( "Thermal contact conductance overflow" );
    matrix[a*m+a] += g; matrix[b*m+b] += g;
    matrix[a*m+b] -= g; matrix[b*m+a] -= g;
  }
  // Symmetric Dirichlet elimination preserves the SPD system.
  for( std::size_t i = 0; i < m; ++i ) if( fixed[i] )
  {
    for( std::size_t j = 0; j < m; ++j ) if( j != i )
    {
      rhs[j] -= matrix[j*m+i]*fixedTemperature[i];
      matrix[j*m+i] = 0.0; matrix[i*m+j] = 0.0;
    }
    matrix[i*m+i] = capacity[i]; rhs[i] = capacity[i]*fixedTemperature[i];
  }
  // Diagonal equilibration considerably improves unequal-capacity cases.
  std::vector< double > scale( m );
  for( std::size_t i = 0; i < m; ++i ) scale[i] = std::sqrt( matrix[i*m+i] );
  for( std::size_t i = 0; i < m; ++i )
  {
    rhs[i] /= scale[i];
    for( std::size_t j = 0; j < m; ++j ) matrix[i*m+j] /= scale[i]*scale[j];
  }
  auto temperature = solveSPD( std::move( matrix ), std::move( rhs ) );
  for( std::size_t i = 0; i < m; ++i ) temperature[i] /= scale[i];
  std::vector< double > groupEnergy( m, 0.0 ), result( n, 0.0 ), freeChange( m, 0.0 ), fixedCapacity( m, 0.0 );
  for( std::size_t i = 0; i < n; ++i ) groupEnergy[index[i]] += nodes[i].bulkEnergy;
  for( auto const & edge : edges ) if( !edge.perfect )
  {
    auto const a = index[edge.a], b = index[edge.b];
    if( a == b ) continue;
    double const heat = dt*edge.conductance*( temperature[b]-temperature[a] );
    groupEnergy[a] += heat; groupEnergy[b] -= heat;
  }
  for( std::size_t i = 0; i < n; ++i )
  {
    auto const r = index[i];
    if( nodes[i].prescribed ) fixedCapacity[r] += nodes[i].capacity;
    else
    {
      result[i] = nodes[i].capacity*( temperature[r]-nodes[i].temperature );
      freeChange[r] += result[i];
    }
  }
  for( std::size_t i = 0; i < n; ++i )
  {
    auto const r = index[i];
    if( fixedCapacity[r] > 0.0 )
    {
      if( nodes[i].prescribed ) result[i] = ( groupEnergy[r]-freeChange[r] )*nodes[i].capacity/fixedCapacity[r];
    }
    else result[i] += ( groupEnergy[r]-freeChange[r] )*nodes[i].capacity/capacity[r];
  }
  return result;
}

inline double basisValue( double const coordinate, bool const cubic )
{
  double const r = std::abs( coordinate );
  if( !cubic ) return std::max( 0.0, 1.0-r );
  return r < 1.0 ? 2.0/3.0-r*r+0.5*r*r*r : ( r < 2.0 ? std::pow( 2.0-r, 3 )/6.0 : 0.0 );
}

inline double basisIntegral( double const lower, double const upper, double const h, bool const cubic )
{
  auto primitive = [=]( double const x )
  {
    double const s = std::abs( x )/h;
    double f;
    if( cubic ) f = s < 1.0 ? 2.0*s/3.0-s*s*s/3.0+s*s*s*s/8.0 :
                               ( s < 2.0 ? 0.5-std::pow( 2.0-s, 4 )/24.0 : 0.5 );
    else { double const r = std::min( 1.0, s ); f = r-0.5*r*r; }
    return std::copysign( h*f, x );
  };
  return upper > lower ? primitive( upper )-primitive( lower ) : 0.0;
}

/** Integral of the linear/cubic grid basis over a plane clipped to its support
 * and physical domain. Grid-aligned planes are exact. Other planes use midpoint
 * surface quadrature; one unique area transfers heat on both sides.
 */
inline double planeNodalArea( Vector normal, Vector const & point, Vector const & h,
                              Vector const & lower, Vector const & upper, int const order,
                              double const cubicFractionA = 0.0, double const cubicFractionB = 0.0 )
{
  if( !normalize( normal ) || order < 2 ) return 0.0;
  double const fraction = std::min( 1.0, std::max( 0.0, 0.5*( cubicFractionA+cubicFractionB ) ) );
  for( int axisIndex = 0; axisIndex < 3; ++axisIndex ) if( std::abs( normal[axisIndex] ) > 1.0-1.0e-12 )
  {
    double const coordinate = dot( normal, point )/normal[axisIndex];
    if( coordinate < lower[axisIndex] || coordinate > upper[axisIndex] ) return 0.0;
    double linear = basisValue( coordinate/h[axisIndex], false );
    double cubic = basisValue( coordinate/h[axisIndex], true );
    for( int k = 0; k < 3; ++k ) if( k != axisIndex )
    {
      linear *= basisIntegral( lower[k], upper[k], h[k], false );
      cubic *= basisIntegral( lower[k], upper[k], h[k], true );
    }
    return ( 1.0-fraction )*linear+fraction*cubic;
  }
  Vector axis = {{ 1.0, 0.0, 0.0 }};
  if( std::abs( normal[0] ) > 0.8 ) axis = {{ 0.0, 1.0, 0.0 }};
  Vector u = cross( normal, axis ); normalize( u );
  Vector const v = cross( normal, u );
  double const supportScale = cubicFractionA+cubicFractionB > 0.0 ? 2.0 : 1.0;
  double const radius = supportScale*std::sqrt( dot( h, h ) );
  double const increment = 2.0*radius/order;
  // Closest point to the node avoids a tangentially shifted integration square.
  double const distance = dot( normal, point );
  double area = 0.0;
  for( int i = 0; i < order; ++i ) for( int j = 0; j < order; ++j )
  {
    double const a = -radius+( i+0.5 )*increment;
    double const b = -radius+( j+0.5 )*increment;
    double linear = 1.0, cubic = 1.0;
    for( int k = 0; k < 3; ++k )
    {
      double const x = distance*normal[k]+a*u[k]+b*v[k];
      if( x < lower[k] || x > upper[k] ) { linear = 0.0; cubic = 0.0; break; }
      double const r = std::abs( x )/h[k];
      linear *= basisValue( r, false );
      cubic *= basisValue( r, true );
    }
    double const weight = ( 1.0-fraction )*linear+fraction*cubic;
    area += weight*increment*increment;
  }
  return area;
}

} // namespace mpmThermal
} // namespace geos
#endif

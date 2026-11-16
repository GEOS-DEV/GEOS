/* SPDX-License-Identifier: LGPL-2.1-only */
#ifndef GEOS_CONSTITUTIVE_THERMOMECHANICALLAW_HPP_
#define GEOS_CONSTITUTIVE_THERMOMECHANICALLAW_HPP_

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace geos
{
namespace constitutive
{
namespace thermomechanics
{
using Vector = std::array< double, 6 >;
using Matrix = std::array< Vector, 6 >;

inline double dot( Vector const & a, Vector const & b )
{
  double r = 0.0;
  for( unsigned i = 0; i < 6; ++i ) r += a[i]*b[i];
  return r;
}

inline Vector matvec( Matrix const & a, Vector const & x )
{
  Vector y{};
  for( unsigned i = 0; i < 6; ++i ) y[i] = dot( a[i], x );
  return y;
}

/** Engineering-strain Voigt convention: stress dot strain is work, without
 * extra factors on shear. Diagonal equilibration avoids a unit-sized pivot test.
 */
inline Matrix compliance( Matrix const & stiffness )
{
  Matrix a = stiffness, result{};
  Vector scale{};
  for( unsigned i = 0; i < 6; ++i )
  {
    if( !std::isfinite( a[i][i] ) || !( a[i][i] > 0.0 ) )
      throw std::domain_error( "Thermomechanical stiffness must be positive definite" );
    scale[i] = std::sqrt( a[i][i] );
  }
  Matrix l{};
  for( unsigned i = 0; i < 6; ++i )
    for( unsigned j = 0; j <= i; ++j )
    {
      double v = a[i][j]/scale[i]/scale[j];
      if( !std::isfinite( v ) || std::abs( a[i][j]/scale[i]/scale[j]-a[j][i]/scale[j]/scale[i] ) > 1.0e-10 )
        throw std::domain_error( "Thermomechanical stiffness must be symmetric" );
      for( unsigned k = 0; k < j; ++k ) v -= l[i][k]*l[j][k];
      if( i == j )
      {
        if( !( v > 1.0e-14 ) ) throw std::domain_error( "Thermomechanical stiffness is singular or indefinite" );
        l[i][j] = std::sqrt( v );
      }
      else l[i][j] = v/l[j][j];
    }
  for( unsigned column = 0; column < 6; ++column )
  {
    Vector x{};
    for( unsigned i = 0; i < 6; ++i )
    {
      double v = i == column ? 1.0/scale[column] : 0.0;
      for( unsigned j = 0; j < i; ++j ) v -= l[i][j]*x[j];
      x[i] = v/l[i][i];
    }
    for( int i = 5; i >= 0; --i )
    {
      double v = x[i];
      for( unsigned j = i+1; j < 6; ++j ) v -= l[j][i]*x[j];
      x[i] = v/l[i][i];
    }
    for( unsigned i = 0; i < 6; ++i ) result[i][column] = x[i]/scale[i];
  }
  return result;
}

/** Constant elastic stiffness, corotational small elastic strains. Density is
 * reference density; all energies are specific (energy/mass). The thermoelastic
 * term follows u = psi + T*s for psi_el = 1/2 eps_e:C:eps_e,
 * eps_e = eps - eps_p - alpha*(T-Tr)*I. It must not be mistaken for heat.
 */
struct EnergyState
{
  double elastic = 0.0;
  double thermoelasticCoefficient = 0.0;
  double retained = 0.0;
  double stored( double const temperature ) const
  { return elastic+thermoelasticCoefficient*temperature+retained; }
};

inline EnergyState state( Matrix const & s, Vector const & stress,
                          double const referenceDensity, double const alpha,
                          double const retained )
{
  if( !( referenceDensity > 0.0 ) || !std::isfinite( referenceDensity ) )
    throw std::domain_error( "Invalid thermomechanical reference density" );
  return { 0.5*dot( stress, matvec( s, stress ) )/referenceDensity,
           alpha*( stress[0]+stress[1]+stress[2] )/referenceDensity, retained };
}

struct WorkIncrement
{
  double total;
  double plastic;
};

/** Midpoint stress work and plastic work using the actual returned stress.
 * For constant stiffness, eps_p = eps_mech - S:delta(sigma). This supports
 * models such as PerfectlyPlastic that do not store a plastic-strain tensor.
 * Old/new stresses and strain MUST use the same corotated frame.
 */
inline WorkIncrement work( Matrix const & s, Vector const & oldStress,
                           Vector const & newStress, Vector const & strain,
                           double const referenceDensity, double const thermalStrain )
{
  Vector mean{}, delta{}, plastic{};
  for( unsigned i = 0; i < 6; ++i )
  {
    mean[i] = 0.5*( oldStress[i]+newStress[i] );
    delta[i] = newStress[i]-oldStress[i];
  }
  Vector const elastic = matvec( s, delta );
  for( unsigned i = 0; i < 6; ++i )
    plastic[i] = strain[i]-elastic[i]-( i < 3 ? thermalStrain : 0.0 );
  double const total = dot( mean, strain )/referenceDensity;
  double const dissipation = dot( mean, plastic )/referenceDensity;
  double const tolerance = 1.0e-10*std::max( { 1.0, std::abs( total ),
    std::abs( dot( mean, elastic )/referenceDensity ) } );
  if( !std::isfinite( total ) || !std::isfinite( dissipation ) || dissipation < -tolerance )
    throw std::domain_error( "Negative plastic work: material is incompatible with the constant-stiffness energy contract" );
  return { total, std::max( 0.0, dissipation ) };
}

/** Safeguarded local energy solve. Every evaluation must replay a constitutive
 * trial from the SAME committed history. Caller explicitly replays the returned
 * root and commits once. A failed solve must not be accepted as a converged step.
 */
template< typename RESIDUAL >
double solveTemperature( double const initial, double const energyTolerance,
                         int const maxIterations, RESIDUAL && residual )
{
  if( !std::isfinite( initial ) || initial < 0.0 || !( energyTolerance > 0.0 ) || maxIterations < 1 )
    throw std::domain_error( "Invalid local thermomechanical solve controls" );
  auto evaluate = [&]( double t )
  {
    double const r = residual( t );
    if( !std::isfinite( r ) ) throw std::domain_error( "Nonfinite thermomechanical energy residual" );
    return r;
  };
  double const f0 = evaluate( initial );
  if( std::abs( f0 ) <= energyTolerance ) return initial;
  double a = initial, b = initial, fa = f0, fb = f0;
  double width = std::max( 1.0, 0.01*initial );
  bool bracketed = false;
  for( int i = 0; i < maxIterations; ++i )
  {
    // Search in the direction predicted by positive heat capacity.
    if( f0 > 0.0 ) { a = std::max( 0.0, initial-width ); fa = evaluate( a ); }
    else { b = initial+width; fb = evaluate( b ); }
    if( std::abs( fa ) <= energyTolerance ) return a;
    if( std::abs( fb ) <= energyTolerance ) return b;
    if( std::signbit( fa ) != std::signbit( fb ) ) { bracketed = true; break; }
    if( f0 > 0.0 && !( a > 0.0 ) ) break;
    width *= 2.0;
  }
  if( !bracketed ) throw std::runtime_error( "Cannot bracket an admissible thermomechanical temperature; reduce the timestep" );
  for( int i = 0; i < maxIterations; ++i )
  {
    double t = ( a*fb-b*fa )/( fb-fa );
    double const margin = 0.1*( b-a );
    if( !std::isfinite( t ) || t <= a+margin || t >= b-margin ) t = 0.5*( a+b );
    double const f = evaluate( t );
    if( std::abs( f ) <= energyTolerance ) return t;
    if( std::signbit( f ) == std::signbit( fa ) ) { a = t; fa = f; }
    else { b = t; fb = f; }
  }
  throw std::runtime_error( "Local thermomechanical energy solve failed; reduce the timestep" );
}
} // namespace thermomechanics
} // namespace constitutive
} // namespace geos
#endif

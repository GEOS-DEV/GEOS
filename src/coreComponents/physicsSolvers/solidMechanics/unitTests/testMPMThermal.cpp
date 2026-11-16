/* SPDX-License-Identifier: LGPL-2.1-only */
// Also runnable without GEOS/TPLs:
// g++ -std=c++17 -I src/coreComponents testMPMThermal.cpp -o testMPMThermal
#include "physicsSolvers/solidMechanics/MPMThermalNumerics.hpp"
#include "constitutive/solid/ThermalLaw.hpp"
#include <iostream>
#include <random>
#include <string>

namespace
{
using namespace geos::mpmThermal;
int checks = 0;
void require( bool const ok, std::string const & message )
{
  ++checks;
  if( !ok ) throw std::runtime_error( message );
}
void near( double const a, double const b, double const tolerance, std::string const & message )
{ require( std::abs( a-b ) <= tolerance, message ); }

void testCaloricPolicies()
{
  geos::constitutive::ConstantThermalLaw law;
  law.capacity = 500.0; law.referenceTemperature = 300.0;
  near( law.energy( 250.0 ), -25000.0, 1.0e-12, "Reference energy may be negative" );
  near( law.temperature( law.energy( 425.0 ) ), 425.0, 1.0e-12, "Constant caloric inversion" );
  geos::constitutive::LinearThermalLaw linear;
  linear.capacity = 500.0; linear.capacitySlope = 0.75; linear.referenceTemperature = 300.0;
  for( double t : { 100.0, 299.0, 300.0, 500.0, 1000.0 } )
  {
    near( linear.temperature( linear.energy( t ) ), t, 1.0e-10, "Nonlinear caloric inversion" );
    double const step = 1.0e-3;
    near( ( linear.energy( t+step )-linear.energy( t-step ) )/( 2*step ), linear.specificHeat( t ),
          1.0e-6, "Heat capacity must be the derivative of specific energy" );
  }
  // A local source changes energy directly and is never diluted by a P2G pass.
  double const mass = 0.002, power = -3.0, duration = 4.0;
  double const u = law.energy( 400.0 )+power*duration/mass;
  near( law.temperature( u ), 388.0, 1.0e-12, "Analytical particle heat loss" );
}

void testBinaryContact()
{
  std::vector< ThermalNode > n( 2 );
  n[0].capacity = 2.0; n[0].temperature = 400.0;
  n[1].capacity = 5.0; n[1].temperature = 280.0;
  double const dt = 0.3, conductance = 7.0;
  auto const e = exchangeGraph( dt, n, { { 0, 1, conductance, false } } );
  double const exact = dt*conductance*( 400.0-280.0 )/( 1.0+dt*conductance*( 0.5+0.2 ) );
  near( e[0], -exact, 1.0e-11, "Backward Euler binary heat transfer" );
  near( e[0]+e[1], 0.0, 1.0e-12, "Binary contact conservation" );
  auto const isolated = exchangeGraph( dt, n, {} );
  near( isolated[0], 0.0, 1.0e-12, "Insulated interface" );
  auto const perfect = exchangeGraph( dt, n, { { 0, 1, 0.0, true } } );
  double const equilibrium = ( 2.0*400.0+5.0*280.0 )/7.0;
  near( n[0].temperature+perfect[0]/n[0].capacity, equilibrium, 1.0e-11, "Perfect contact temperature A" );
  near( n[1].temperature+perfect[1]/n[1].capacity, equilibrium, 1.0e-11, "Perfect contact temperature B" );
}

void testMultifieldAndReservoir()
{
  std::vector< ThermalNode > n( 3 );
  for( auto & node : n ) node.capacity = 1.0;
  n[0].temperature = 500.0; n[1].temperature = 350.0; n[2].temperature = 200.0;
  auto const excluded = exchangeGraph( 2.0, n, { { 0, 2, 1.0, false } } );
  near( excluded[1], 0.0, 1.0e-12, "An excluded A-B interface must transfer no heat" );
  auto const e = exchangeGraph( 2.0, n, { { 0, 1, 1.0, false }, { 1, 2, 2.0, false }, { 0, 2, 0.5, false } } );
  near( e[0]+e[1]+e[2], 0.0, 1.0e-11, "Three-field graph conservation" );
  for( std::size_t i = 0; i < n.size(); ++i ) require( n[i].temperature+e[i] >= 200.0 && n[i].temperature+e[i] <= 500.0,
                                                    "Implicit graph may not overshoot initial temperatures" );
  n.resize( 2 ); n[0].temperature = 400.0; n[0].prescribed = true; n[0].prescribedTemperature = 400.0;
  n[1].temperature = 200.0;
  auto const bath = exchangeGraph( 2.0, n, { { 0, 1, 3.0, false } } );
  double const final = ( 200.0+6.0*400.0 )/7.0;
  near( n[1].temperature+bath[1], final, 1.0e-10, "Isothermal reservoir exchange" );
  near( bath[0]+bath[1], 0.0, 1.0e-12, "Reservoir physical exchange is conservative before thermostat reaction" );
  // Two prescribed particles can absorb bulk/source heat without changing T.
  n[0].bulkEnergy = 20.0; n[1].bulkEnergy = -20.0;
  auto const withBulk = exchangeGraph( 0.1, n, { { 0, 1, 3.0, false } } );
  near( withBulk[0]+withBulk[1], 0.0, 1.0e-11, "Bulk and contact energy accounting" );
}

void testCapacityWeightedTransfer()
{
  // Unequal masses and heat capacities, two particles sharing two grid nodes.
  double const mass[2] = { 2.0, 5.0 }, c[2] = { 3.0, 7.0 };
  double const shape[2][2] = { { 0.8, 0.2 }, { 0.3, 0.7 } };
  double capacity[2] = {};
  for( int p = 0; p < 2; ++p ) for( int i = 0; i < 2; ++i ) capacity[i] += mass[p]*c[p]*shape[p][i];
  double const energy[2] = { -12.0, 12.0 };
  double particleTotal = 0.0;
  for( int p = 0; p < 2; ++p ) for( int i = 0; i < 2; ++i )
    particleTotal += mass[p]*gatherSpecificEnergy( c[p], shape[p][i], energy[i], capacity[i] );
  near( particleTotal, 0.0, 1.0e-12, "P2G/G2P must conserve extensive energy with unequal specific heats" );
  // Per-original-particle power must survive subdivision.
  double const originalMass = 8.0, originalPower = 16.0;
  double totalPower = 0.0;
  for( double childMass : { 1.0, 2.0, 5.0 } ) totalPower += originalPower*childMass/originalMass;
  near( totalPower, originalPower, 1.0e-12, "Subdivision must conserve source power" );
}

void testAreas()
{
  Vector const h = {{ 1.0, 1.0, 1.0 }}, lower = {{ -1.0, -1.0, -1.0 }}, upper = {{ 1.0, 1.0, 1.0 }};
  double const full = planeNodalArea( {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, h, lower, upper, 240 );
  near( full, 1.0, 5.0e-4, "Integrated planar nodal area" );
  Vector const halfUpper = {{ 1.0, 1.0, 0.0 }};
  double const half = planeNodalArea( {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, h, lower, halfUpper, 240 );
  near( half, 0.5, 3.0e-4, "Physical-domain clipping prevents doubled plane-strain thickness" );
  double const offset = planeNodalArea( {{ 1.0, 0.0, 0.0 }}, {{ 0.25, 0.0, 0.0 }}, h, lower, upper, 240 );
  near( offset, 0.75, 5.0e-4, "Grid-offset interface area" );
  Vector const lo2 = {{ -2.0, -2.0, -2.0 }}, hi2 = {{ 2.0, 2.0, 2.0 }};
  double const cubic = planeNodalArea( {{ 1.0, 0.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, h, lo2, hi2, 240, 1.0, 1.0 );
  near( cubic, 2.0/3.0, 5.0e-4, "Cubic B-spline interface basis" );
  double const diagonal = planeNodalArea( {{ 1.0, 1.0, 0.0 }}, {{ 0.0, 0.0, 0.0 }}, h, lower, upper, 240 );
  near( diagonal, 2.0*std::sqrt( 2.0 )/3.0, 5.0e-4, "Oblique planar nodal area" );
}

void testBulkConduction()
{
  std::vector< Vector > const gradient = { {{ -2.0, 0.0, 0.0 }}, {{ 2.0, 0.0, 0.0 }} };
  auto const hot = assembleFourierParticle( 0.5, 3.0, gradient, { 400.0, 300.0 } );
  near( hot.power[0], -600.0, 1.0e-12, "Fourier weak-form power on a hot node" );
  near( hot.power[0]+hot.power[1], 0.0, 1.0e-12, "Bulk conduction conserves energy" );
  auto const uniform = assembleFourierParticle( 0.5, 3.0, gradient, { 350.0, 350.0 } );
  near( uniform.power[0], 0.0, 1.0e-12, "Uniform temperature is preserved" );
  auto const adiabatic = assembleFourierParticle( 0.5, 0.0, gradient, { 400.0, 300.0 } );
  near( adiabatic.power[0], 0.0, 1.0e-12, "Zero conductivity has no grid-remap diffusion" );
  // Each quadrature contribution can be assembled on a different MPI rank.
  auto const rank0 = assembleFourierParticle( 0.2, 3.0, gradient, { 400.0, 300.0 } );
  auto const rank1 = assembleFourierParticle( 0.3, 3.0, gradient, { 400.0, 300.0 } );
  near( rank0.power[0]+rank1.power[0], hot.power[0], 1.0e-12, "Partitioned assembly equals serial assembly" );
  // A two-node lumped bar has d(Ta-Tb)/dt=-2*G/C*(Ta-Tb).
  double ta = 400.0, tb = 300.0;
  double const dt = 1.0e-4, end = 0.01, capacity = 2.0;
  for( int step = 0; step < 100; ++step )
  {
    auto const rhs = assembleFourierParticle( 0.5, 3.0, gradient, { ta, tb } );
    ta += dt*rhs.power[0]/capacity; tb += dt*rhs.power[1]/capacity;
  }
  near( ta-tb, 100*std::exp( -6.0*end ), 0.002, "Bulk transient approaches the analytical exponential" );
  near( ta+tb, 700.0, 1.0e-10, "Transient bar energy balance" );
}

void testSlabResistance()
{
  // Stationary unit-length bar, k=rho*c=1, unit area, two thermal fields,
  // R''=0.2, endpoint reservoirs at 400/200 K. Assemble the complete operator
  // cycle using the same Fourier, contact, and energy-gather routines as GEOS.
  int const cells = 16;
  double const h = 1.0/cells, dt = 0.1*h*h, resistance = 0.2;
  std::vector< double > particleTemperature( cells, 300.0 );
  particleTemperature.front() = 400.0; particleTemperature.back() = 200.0;
  double hotBathPower = 0.0;
  for( int step = 0; step < static_cast< int >( 5.0/dt ); ++step )
  {
    std::vector< std::array< double, 2 > > capacity( cells+1 ), moment( cells+1 ),
      temperature( cells+1 ), power( cells+1 ), energy( cells+1 );
    for( int p = 0; p < cells; ++p )
    {
      int const field = p < cells/2 ? 0 : 1;
      for( int j = 0; j < 2; ++j )
      {
        capacity[p+j][field] += 0.5*h;
        moment[p+j][field] += 0.5*h*particleTemperature[p];
      }
    }
    for( int i = 0; i <= cells; ++i ) for( int f = 0; f < 2; ++f )
      if( capacity[i][f] > 0.0 ) temperature[i][f] = moment[i][f]/capacity[i][f];
    for( int p = 0; p < cells; ++p )
    {
      int const f = p < cells/2 ? 0 : 1;
      auto const rhs = assembleFourierParticle( h, 1.0, { {{ -1/h, 0.0, 0.0 }}, {{ 1/h, 0.0, 0.0 }} },
                                                { temperature[p][f], temperature[p+1][f] } );
      for( int j = 0; j < 2; ++j ) power[p+j][f] += rhs.power[j];
    }
    for( int i = 0; i <= cells; ++i ) for( int f = 0; f < 2; ++f ) energy[i][f] = dt*power[i][f];
    int const interface = cells/2;
    std::vector< ThermalNode > state( 2 );
    for( int f = 0; f < 2; ++f )
    {
      state[f].capacity = capacity[interface][f]; state[f].temperature = temperature[interface][f];
      state[f].bulkEnergy = energy[interface][f];
    }
    auto const contact = exchangeGraph( dt, state, { { 0, 1, 1.0/resistance, false } } );
    energy[interface][0] = contact[0]; energy[interface][1] = contact[1];
    for( int p = 0; p < cells; ++p )
    {
      int const f = p < cells/2 ? 0 : 1;
      for( int j = 0; j < 2; ++j )
        particleTemperature[p] += gatherSpecificEnergy( 1.0, 0.5, energy[p+j][f], capacity[p+j][f] );
    }
    hotBathPower = h*( 400.0-particleTemperature.front() )/dt;
    particleTemperature.front() = 400.0; particleTemperature.back() = 200.0;
  }
  near( hotBathPower, 200.0/( 1.0+resistance ), 0.02, "Two-material slab analytical series resistance" );
}

void testRandomGraphs()
{
  std::mt19937 generator( 1729 );
  std::uniform_real_distribution< double > uniform( 0.1, 1.0 );
  for( int trial = 0; trial < 100; ++trial )
  {
    std::vector< ThermalNode > n( 6 ); std::vector< ThermalEdge > edges;
    for( auto & node : n ) { node.capacity = 100*uniform( generator ); node.temperature = 1000*uniform( generator ); }
    for( std::size_t a = 0; a < n.size(); ++a ) for( std::size_t b = a+1; b < n.size(); ++b )
      edges.push_back( { a, b, 100*uniform( generator ), false } );
    auto const e = exchangeGraph( 3.0, n, edges );
    double sum = 0.0;
    double lo = n.front().temperature, hi = lo;
    for( auto const & node : n ) { lo = std::min( lo, node.temperature ); hi = std::max( hi, node.temperature ); }
    for( std::size_t i = 0; i < n.size(); ++i )
    {
      sum += e[i];
      double const t = n[i].temperature+e[i]/n[i].capacity;
      require( t >= lo-1.0e-10 && t <= hi+1.0e-10, "Random graph maximum principle" );
    }
    near( sum, 0.0, 1.0e-8, "Random graph conservation" );
  }
  bool rejected = false;
  try { exchangeGraph( -1.0, {}, {} ); } catch( std::invalid_argument const & ) { rejected = true; }
  require( rejected, "Invalid time step must be rejected" );
  auto const tinyDt = exchangeGraph( 1.0e-16, { ThermalNode{ 1.0, 300.0, 0.0, false, 0.0 } }, {} );
  near( tinyDt[0], 0.0, 1.0e-10, "A small positive time step is not zero" );
}
}

int main()
{
  try
  {
    testCaloricPolicies(); testBinaryContact(); testMultifieldAndReservoir();
    testCapacityWeightedTransfer(); testAreas(); testRandomGraphs(); testBulkConduction(); testSlabResistance();
    std::cout << "MPM thermal numerical checks passed: " << checks << '\n';
  }
  catch( std::exception const & error ) { std::cerr << error.what() << '\n'; return 1; }
  return 0;
}

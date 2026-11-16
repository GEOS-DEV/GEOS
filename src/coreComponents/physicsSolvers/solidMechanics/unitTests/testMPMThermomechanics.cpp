/* SPDX-License-Identifier: LGPL-2.1-only */
// Standalone: g++ -std=c++17 -I src/coreComponents testMPMThermomechanics.cpp -o testMPMThermomechanics
#include "constitutive/solid/ThermalLaw.hpp"
#include <iostream>
#include <string>

namespace thermo = geos::constitutive::thermomechanics;
namespace
{
int checks = 0;
void check( bool const condition, char const * message )
{ ++checks; if( !condition ) throw std::runtime_error( message ); }
void near( double a, double b, double tol, char const * message )
{ check( std::abs( a-b ) <= tol, message ); }

thermo::Matrix stiffness()
{
  thermo::Matrix c{};
  double const k = 1000.0, g = 500.0, lambda = k-2.0*g/3.0;
  for( unsigned i = 0; i < 3; ++i )
    for( unsigned j = 0; j < 3; ++j ) c[i][j] = lambda+( i == j ? 2.0*g : 0.0 );
  for( unsigned i = 3; i < 6; ++i ) c[i][i] = g;
  return c;
}

double mises( thermo::Vector const & s )
{
  double const mean = ( s[0]+s[1]+s[2] )/3.0;
  return std::sqrt( 1.5*( ( s[0]-mean )*( s[0]-mean )+( s[1]-mean )*( s[1]-mean )+
    ( s[2]-mean )*( s[2]-mean )+2.0*( s[3]*s[3]+s[4]*s[4]+s[5]*s[5] ) ) );
}

struct Particle
{
  thermo::Vector stress{};
  double t = 300.0, u = 0.0, retained = 0.0, plastic = 0.0;
};

// Independent radial-return material for exercising the production energy and
// root helpers. GEOS constitutive dispatch is tested separately when TPLs exist.
Particle advance( Particle const & old, thermo::Vector const & strain, double beta,
                  double softening, bool lagged = false, bool fixed = false,
                  double alpha = 0.0, double capacitySlope = 0.0 )
{
  auto const c = stiffness(), s = thermo::compliance( c );
  geos::constitutive::LinearThermalLaw law;
  law.capacity = 1.0; law.referenceTemperature = 300.0; law.capacitySlope = capacitySlope;
  Particle trial;
  double work = 0.0;
  auto evaluate = [&]( double t )
  {
    // Each trial starts with OLD history, including a deliberately observable
    // cumulative plastic-work history. Iteration count cannot multiply heating.
    trial = old;
    trial.t = t;
    auto mechanical = strain;
    for( unsigned i = 0; i < 3; ++i ) mechanical[i] -= alpha*( t-old.t );
    auto ds = thermo::matvec( c, mechanical );
    for( unsigned i = 0; i < 6; ++i ) trial.stress[i] += ds[i];
    double const q = mises( trial.stress ), yield = 10.0*std::exp( -softening*( t-300.0 ) );
    if( q > yield )
    {
      double const mean = ( trial.stress[0]+trial.stress[1]+trial.stress[2] )/3.0;
      for( unsigned i = 0; i < 6; ++i )
        trial.stress[i] = ( i < 3 ? mean : 0.0 )+( trial.stress[i]-( i < 3 ? mean : 0.0 ) )*yield/q;
    }
    auto const w = thermo::work( s, old.stress, trial.stress, strain, 1.0, alpha*( t-old.t ) );
    trial.retained = old.retained+( 1.0-beta )*w.plastic;
    trial.plastic = old.plastic+w.plastic;
    auto const state = thermo::state( s, trial.stress, 1.0, alpha, trial.retained );
    trial.u = law.energy( t )+state.stored( t );
    work = w.total;
    return trial.u-old.u-work;
  };
  if( fixed )
  {
    evaluate( old.t );
    // Bath is evaluated independently by the test from this accepted state.
  }
  else if( lagged )
  {
    evaluate( old.t );
    auto const state = thermo::state( s, trial.stress, 1.0, alpha, trial.retained );
    double const target = old.u+work;
    trial.t = thermo::solveTemperature( old.t, 1.0e-11, 100,
      [&]( double t ) { return law.energy( t )+state.stored( t )-target; } );
    trial.u = target;
  }
  else
  {
    double const t = thermo::solveTemperature( old.t, 1.0e-11, 100, evaluate );
    evaluate( t );
    near( trial.u-old.u, work, 2.0e-10, "First law: total work counted once" );
  }
  return trial;
}

class TestCaloricMaterial : public geos::constitutive::ThermalMaterial
{
public:
  geos::constitutive::LinearThermalLaw law;
  double thermalSpecificHeat( double t ) const override { return law.specificHeat( t ); }
  double thermalConductivity( double t ) const override { return law.thermalConductivity( t ); }
  double thermalEnergy( double t ) const override { return law.energy( t ); }
  double thermalTemperature( double u ) const override { return law.temperature( u ); }
  bool solverThermalExpansion() const override { return false; }
  bool supportsThermomechanicalEnergy() const override { return true; }
  thermo::Matrix thermomechanicalCompliance( std::size_t ) const override { return thermo::compliance( stiffness() ); }
  double plasticHeatFraction() const override { return 0.9; }
  void prepareThermomechanicalTrial( std::size_t, double ) override {}
};

void testTotalEnergyInversion()
{
  TestCaloricMaterial material;
  material.law.referenceTemperature = 300.0;
  material.law.capacity = 500.0;
  material.law.capacitySlope = 0.1;
  thermo::EnergyState state{ 200000.0, -0.01, 40000.0 };
  for( double t : { 100.0, 300.0, 600.0 } )
  {
    double const u = material.internalEnergy( t, state );
    near( material.temperatureFromInternalEnergy( u, 300.0, state ), t, 1.0e-7,
          "Total-energy inversion excludes stored elastic and retained energy" );
    near( material.stateHeatCapacity( t, state ), material.thermalSpecificHeat( t )-0.01, 1.0e-10,
          "Transport capacity differentiates frozen-state total energy" );
  }
  thermo::EnergyState none{};
  double const negative = material.internalEnergy( 100.0, none );
  check( negative < 0.0, "Reference-relative total energy can be negative" );
  near( material.temperatureFromInternalEnergy( negative, 300.0, none ), 100.0, 1.0e-7,
        "Negative total energy is not clamped to zero" );
}

void testElasticAndWork()
{
  auto const c = stiffness(), s = thermo::compliance( c );
  for( unsigned i = 0; i < 6; ++i )
    for( unsigned j = 0; j < 6; ++j )
    {
      double product = 0.0;
      for( unsigned k = 0; k < 6; ++k ) product += c[i][k]*s[k][j];
      near( product, i == j ? 1.0 : 0.0, 1.0e-12, "Elastic compliance inversion" );
    }
  Particle p;
  thermo::Vector strain{}; strain[0] = 0.001;
  auto loaded = advance( p, strain, 0.9, 0.0 );
  near( loaded.t, 300.0, 1.0e-10, "Elastic loading is not plastic heat" );
  near( loaded.plastic, 0.0, 1.0e-12, "Elastic loading has no plastic work" );
  strain[0] = -0.001;
  auto unloaded = advance( loaded, strain, 0.9, 0.0 );
  near( unloaded.u, 0.0, 1.0e-10, "Elastic energy recovered on unload" );
  near( unloaded.t, 300.0, 1.0e-10, "Elastic unload does not create heat" );
  thermo::Vector shear{}; shear[5] = 0.002;
  auto sheared = advance( p, shear, 0.9, 0.0 );
  near( sheared.u, 0.5*500.0*0.002*0.002, 1.0e-12, "Engineering shear work has correct factor" );
  thermo::Matrix tiny = c;
  for( auto & row : tiny ) for( double & v : row ) v *= 1.0e-20;
  check( std::isfinite( thermo::compliance( tiny )[0][0] ), "Unit-independent stiffness pivots" );
}

void testPlasticCoupling()
{
  thermo::Vector strain{}; strain[5] = 0.2;
  Particle p;
  for( double beta : { 0.0, 0.5, 0.9, 1.0 } )
  {
    auto q = advance( p, strain, beta, 0.05 );
    check( q.plastic > 0.0, "Plastic work is recorded" );
    near( q.retained, ( 1.0-beta )*q.plastic, 1.0e-11, "Unconverted plastic work stays in internal energy" );
    near( q.t-300.0, beta*q.plastic, 2.0e-10, "Plastic work heats caloric component once" );
    near( mises( q.stress ), 10.0*std::exp( -0.05*( q.t-300.0 ) ), 1.0e-10, "Return lies on final-temperature yield surface" );
  }
  auto coupled = advance( p, strain, 0.9, 0.5 );
  auto lagged = advance( p, strain, 0.9, 0.5, true );
  check( mises( coupled.stress ) < mises( lagged.stress ), "Lagged yield surface misses same-step softening" );
  double priorError = 1.0e99;
  for( unsigned n : { 10u, 40u, 160u, 640u } )
  {
    Particle a, b;
    thermo::Vector ds{}; ds[5] = strain[5]/n;
    for( unsigned i = 0; i < n; ++i )
    {
      a = advance( a, ds, 0.9, 0.05 );
      b = advance( b, ds, 0.9, 0.05, true );
    }
    double const error = std::abs( mises( a.stress )-mises( b.stress ) );
    check( error < priorError, "Lagged/coupled difference decreases under timestep refinement" );
    priorError = error;
  }
  auto nonlinear = advance( p, strain, 0.9, 0.05, false, false, 0.0, 0.5 );
  near( ( nonlinear.t-300.0 )+0.25*std::pow( nonlinear.t-300.0, 2 ), 0.9*nonlinear.plastic,
        2.0e-10, "Nonlinear caloric law participates in local solve" );
  auto fixed = advance( p, strain, 0.9, 0.05, false, true );
  auto const w = thermo::work( thermo::compliance( stiffness() ), p.stress, fixed.stress, strain, 1.0, 0.0 );
  double const reservoir = fixed.u-p.u-w.total;
  near( reservoir, -0.9*fixed.plastic, 1.0e-10, "Isothermal reservoir removes plastic heat" );
}

void testThermoelasticityAndGuards()
{
  Particle p;
  thermo::Vector compression{};
  compression[0] = compression[1] = compression[2] = -1.0e-4;
  auto q = advance( p, compression, 0.9, 0.0, false, false, 1.0e-5 );
  check( q.t > p.t, "Adiabatic elastic compression causes reversible thermoelastic heating" );
  near( q.plastic, 0.0, 1.0e-10, "Thermoelastic heating is not plastic dissipation" );
  bool caught = false;
  try { thermo::solveTemperature( 300.0, 1.0e-12, 2, []( double ) { return 1.0; } ); }
  catch( std::runtime_error const & ) { caught = true; }
  check( caught, "Unbracketed energy solve cannot silently commit" );
  caught = false;
  try { thermo::solveTemperature( 300.0, 1.0e-12, 2, []( double ) { return std::numeric_limits< double >::quiet_NaN(); } ); }
  catch( std::domain_error const & ) { caught = true; }
  check( caught, "Nonfinite energy residual rejected" );
  thermo::Matrix bad = stiffness(); bad[0][0] = -1.0;
  caught = false;
  try { thermo::compliance( bad ); }
  catch( std::domain_error const & ) { caught = true; }
  check( caught, "Indefinite stored-energy law rejected" );
  int history = 0, evaluations = 0;
  auto root = thermo::solveTemperature( 300.0, 1.0e-12, 100, [&]( double t )
  { history = 0; ++history; ++evaluations; return t-340.0; } );
  near( root, 340.0, 1.0e-10, "Safeguarded energy root" );
  check( history == 1 && evaluations > 1, "Transactional replay does not multiply constitutive history" );
}
}

int main()
{
  try
  {
    testTotalEnergyInversion(); testElasticAndWork(); testPlasticCoupling(); testThermoelasticityAndGuards();
    std::cout << checks << " thermomechanical checks passed\n";
    return 0;
  }
  catch( std::exception const & e ) { std::cerr << e.what() << '\n'; return 1; }
}

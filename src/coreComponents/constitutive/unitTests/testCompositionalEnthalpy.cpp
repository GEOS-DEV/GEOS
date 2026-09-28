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

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/models/CompositionalEnthalpy.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/EquationOfState.hpp"
#include "constitutive/fluid/multifluid/compositional/parameters/HeatCapacityCoefficients.hpp"

#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

template< int NC >
using EnthalpyData = std::tuple<
  real64 const,       // pressure
  real64 const,       // temperature
  Feed< NC > const,   // phase composition
  integer const,      // phase
  real64 const,       // expected enthalpy
  real64 const        // expected Joule-Thomson coefficient (K/Pa)
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static constexpr integer numPhases = 2;
  static constexpr integer numComps = 4;
  static constexpr integer numCoeffs = 5;

  static std::unique_ptr< TestFluid< numComps > > createFluid()
  {
    return TestFluid< numComps >::create( {Fluid::CO2, Fluid::H2, Fluid::CH4, Fluid::C8H18} );
  }

  static void populateCoefficients( HeatCapacityCoefficients * coefficients )
  {
    coefficients->m_referenceEnthalpy.resize( numPhases, numComps );
    coefficients->m_referenceEnthalpy.zero();

    coefficients->m_coefficients.resize( numComps, numCoeffs );
    std::array< real64, numComps * numCoeffs > coefficientsData {
      3.45314594e+01, 7.07764122e-02, 1.69716907e-04, -1.40256545e-06, 1.83892934e-09,
      2.87501790e+01, 1.01246159e-02, -5.85965230e-05, 1.30968549e-07, -9.36636594e-11,
      3.55287791e+01, 4.22215766e-02, 1.34958295e-04, -4.13365988e-07, 4.05670019e-10,
      2.01824087e+02, 3.07940012e-01, 6.82371480e-04, -1.07019612e-06, -3.48458322e-11
    };
    for( int ic = 0; ic < numComps; ++ic )
    {
      for( int j = 0; j < 5; ++j )
      {
        coefficients->m_coefficients( ic, j ) = coefficientsData[ic*numCoeffs + j];
      }
    }
    coefficients->m_phaseTypes.resize( numPhases );
    for( PhaseType const phaseType : {PhaseType::LIQUID, PhaseType::VAPOUR} )
    {
      integer const phase = static_cast< integer >(phaseType);
      coefficients->m_phaseTypes[phase] = phase;
    }
  }
};

template<>
struct FluidData< 6 >
{
  static constexpr integer numPhases = 2;
  static constexpr integer numComps = 6;
  static constexpr integer numCoeffs = 5;

  static std::unique_ptr< TestFluid< numComps > > createFluid()
  {
    return TestFluid< numComps >::create( {Fluid::CO2, Fluid::N2, Fluid::H2O, Fluid::CH4, Fluid::H2S, Fluid::C5H12} );
  }

  static void populateCoefficients( HeatCapacityCoefficients * coefficients )
  {
    coefficients->m_referenceEnthalpy.resize( numPhases, numComps );
    coefficients->m_referenceEnthalpy.zero();

    coefficients->m_coefficients.resize( numComps, numCoeffs );
    std::array< real64, numComps * numCoeffs > coefficientsData {
      3.40933377e+01, 7.52742132e-02, 2.00271751e-04, -1.61900936e-06, 2.12199565e-09,
      2.89614125e+01, 2.03166820e-03, 2.47038198e-05, -9.54674193e-08, 1.15955044e-10,
      2.69032982e+01, 5.51993739e-02, 6.31405927e-04, -3.66166336e-06, 4.56291585e-09,
      3.54137115e+01, 4.34392612e-02, 1.42852598e-04, -4.70393972e-07, 4.80614694e-10,
      3.18873193e+01, 3.25765666e-02, 2.20667710e-04, -1.29890488e-06, 1.63371768e-09,
      1.21004526e+02, 2.47982404e-01, 7.27115529e-04, -2.64999262e-06, 2.52476069e-09
    };
    for( int ic = 0; ic < numComps; ++ic )
    {
      for( int j = 0; j < 5; ++j )
      {
        coefficients->m_coefficients( ic, j ) = coefficientsData[ic*numCoeffs + j];
      }
    }
    coefficients->m_phaseTypes.resize( numPhases );
    coefficients->m_phaseTypes[static_cast< integer >(PhaseType::LIQUID)] = PhaseType::LIQUID;
    coefficients->m_phaseTypes[static_cast< integer >(PhaseType::VAPOUR)] = PhaseType::VAPOUR;
  }
};

template< EquationOfStateType EOS_TYPE, int NUM_COMP >
class CompositionalEnthalpyTestFixture : public ::testing::TestWithParam< EnthalpyData< NUM_COMP > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NUM_COMP;
  static constexpr int numDofs = NUM_COMP + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;

public:
  using DataType = EnthalpyData< numComps >;

public:
  CompositionalEnthalpyTestFixture()
    : m_fluid( FluidData< numComps >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();
    m_parameters = CompositionalEnthalpy::createParameters( std::make_unique< ModelParameters >() );

    auto equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EOS_TYPE );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto heatCapacityCoefficients = const_cast< HeatCapacityCoefficients * >(m_parameters->get< HeatCapacityCoefficients >());
    FluidData< numComps >::populateCoefficients( heatCapacityCoefficients );

    string const name = GEOS_FMT( "PhaseEnthalpy{}{}", eosName, numComps );
    m_enthalpy = std::make_unique< CompositionalEnthalpy >( name, componentProperties, 0, *m_parameters );
  }

  ~CompositionalEnthalpyTestFixture() = default;

  void testEnthalpyValues( DataType const & data, bool useMass = false )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, std::get< 2 >( data ));
    setPhase( std::get< 3 >( data ) );
    real64 const & expectedEnthalpy = std::get< 4 >( data );
    real64 const & expectedJouleThomson = std::get< 5 >( data );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = m_enthalpy->createKernelWrapper();

    real64 enthalpy = 0.0;
    stackArray1d< real64, numDofs > enthalpyDerivs( numDofs );

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           enthalpy,
                           enthalpyDerivs.toSlice(),
                           useMass );

    // Approximate JT coefficient
    real64 const jouleThomson = -enthalpyDerivs[Deriv::dP] / enthalpyDerivs[Deriv::dT];

    checkRelativeError( enthalpy, expectedEnthalpy, relTol, absTol );
    checkRelativeError( jouleThomson, expectedJouleThomson, relTol, absTol );
  }

  void testEnthalpyDerivatives( DataType const & data, bool useMass )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< numComps >::createArray( phaseComposition, std::get< 2 >( data ));
    setPhase( std::get< 3 >( data ) );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = m_enthalpy->createKernelWrapper();

    real64 enthalpy = 0.0;
    stackArray2d< real64, 2*numDofs > derivSpace( 2, numDofs );
    arraySlice1d< real64 > enthalpyDerivs = derivSpace[0];
    arraySlice1d< real64 > tempDerivs = derivSpace[1];

    kernelWrapper.compute( componentProperties,
                           pressure,
                           temperature,
                           phaseComposition.toSliceConst(),
                           enthalpy,
                           enthalpyDerivs,
                           useMass );
    // Compare against numerical derivatives
    // -- Pressure derivative
    real64 const dp = 1.0e-4 * pressure;
    internal::testNumericalDerivative(
      pressure, dp, enthalpyDerivs[Deriv::dP],
      [&]( real64 const p ){
      real64 displacedEnthalpy = 0.0;
      kernelWrapper.compute( componentProperties,
                             p,
                             temperature,
                             phaseComposition.toSliceConst(),
                             displacedEnthalpy,
                             tempDerivs,
                             useMass );
      return displacedEnthalpy;
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-5 * temperature;
    internal::testNumericalDerivative(
      temperature, dT, enthalpyDerivs[Deriv::dT],
      [&]( real64 const t ){
      real64 displacedEnthalpy = 0.0;
      kernelWrapper.compute( componentProperties,
                             pressure,
                             t,
                             phaseComposition.toSliceConst(),
                             displacedEnthalpy,
                             tempDerivs,
                             useMass );
      return displacedEnthalpy;
    } );

    // -- Composition derivatives
    real64 constexpr dz = 1.0e-6;
    for( integer ic = 0; ic < numComps; ++ic )
    {
      internal::testNumericalDerivative(
        0.0, dz, enthalpyDerivs[Deriv::dC+ic],
        [&]( real64 const z ){
        real64 const z_old = phaseComposition[ic];
        phaseComposition[ic] += z;
        real64 displacedEnthalpy = 0.0;
        kernelWrapper.compute( componentProperties,
                               pressure,
                               temperature,
                               phaseComposition.toSliceConst(),
                               displacedEnthalpy,
                               tempDerivs,
                               useMass );
        phaseComposition[ic] = z_old;
        return displacedEnthalpy;
      } );
    }
  }

protected:
  void setPhase( integer const phaseIndex )
  {
    ComponentProperties const & componentProperties = m_fluid->getComponentProperties();
    string const name = m_enthalpy->functionName();
    m_enthalpy = std::make_unique< CompositionalEnthalpy >( name, componentProperties, phaseIndex, *m_parameters );
  }

protected:
  std::unique_ptr< CompositionalEnthalpy > m_enthalpy{};
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PengRobinson_4 = CompositionalEnthalpyTestFixture< EquationOfStateType::PengRobinson, 4 >;
using SoaveRedlichKwong_6 = CompositionalEnthalpyTestFixture< EquationOfStateType::SoaveRedlichKwong, 6 >;
using SoreideWhitson_6 = CompositionalEnthalpyTestFixture< EquationOfStateType::SoreideWhitson, 6 >;

TEST_P( PengRobinson_4, testEnthalpyValues )
{
  testEnthalpyValues( GetParam() );
}
TEST_P( SoaveRedlichKwong_6, testEnthalpyValues )
{
  testEnthalpyValues( GetParam() );
}
TEST_P( SoreideWhitson_6, testEnthalpyValues )
{
  testEnthalpyValues( GetParam(), true );
}

TEST_P( PengRobinson_4, testEnthalpyDerivatives )
{
  testEnthalpyDerivatives( GetParam(), false );
  testEnthalpyDerivatives( GetParam(), true );
}
TEST_P( SoaveRedlichKwong_6, testEnthalpyDerivatives )
{
  testEnthalpyDerivatives( GetParam(), false );
  testEnthalpyDerivatives( GetParam(), true );
}
TEST_P( SoreideWhitson_6, testEnthalpyDerivatives )
{
  testEnthalpyDerivatives( GetParam(), false );
  testEnthalpyDerivatives( GetParam(), true );
}

/* UNCRUSTIFY-OFF */

// Test data

constexpr integer LIQ = static_cast<integer>(PhaseType::LIQUID);
constexpr integer VAP = static_cast<integer>(PhaseType::VAPOUR);

INSTANTIATE_TEST_SUITE_P(
  CompositionalEnthalpyTest, PengRobinson_4,
  ::testing::ValuesIn< typename PengRobinson_4::DataType >( {
    {1.01325e+05, 283.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -5.554727e+02,  1.329174e-05},
    {1.00000e+06, 283.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -9.777036e+02,  1.354551e-05},
    {1.00000e+07, 283.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ, -1.258949e+04,  6.097076e-07},
    {1.01325e+05, 298.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -4.172063e+01,  1.189416e-05},
    {1.00000e+06, 298.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -4.271945e+02,  1.203740e-05},
    {1.00000e+07, 298.15, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ, -1.075604e+04,  1.246819e-06},
    {1.01325e+05, 550.00, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP,  1.079339e+04,  2.795143e-06},
    {1.00000e+06, 550.00, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP,  1.067299e+04,  2.740228e-06},
    {1.00000e+07, 550.00, {1.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00}, VAP,  9.544159e+03,  2.175503e-06},
    {1.01325e+05, 283.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -4.305421e+02,  1.699485e-07},
    {1.00000e+06, 283.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -4.346280e+02,  1.478298e-07},
    {1.00000e+07, 283.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ, -4.479027e+02, -3.415949e-08},
    {1.01325e+05, 298.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -3.809923e-01,  1.296216e-07},
    {1.00000e+06, 298.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, VAP, -3.474631e+00,  1.096517e-07},
    {1.00000e+07, 298.15, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ, -9.208397e+00, -5.550701e-08},
    {1.01325e+05, 550.00, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ,  7.363162e+03, -2.039014e-07},
    {1.00000e+06, 550.00, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ,  7.368612e+03, -2.098955e-07},
    {1.00000e+07, 550.00, {0.000000e+00, 1.000000e+00, 0.000000e+00, 0.000000e+00}, LIQ,  7.431415e+03, -2.626841e-07},
    {1.01325e+05, 283.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -5.482831e+02,  5.621645e-06},
    {1.00000e+06, 283.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -7.262985e+02,  5.522977e-06},
    {1.00000e+07, 283.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -2.519542e+03,  3.642260e-06},
    {1.01325e+05, 298.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -1.838859e+01,  5.095639e-06},
    {1.00000e+06, 298.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -1.819651e+02,  4.994247e-06},
    {1.00000e+07, 298.15, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP, -1.781165e+03,  3.367126e-06},
    {1.01325e+05, 550.00, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP,  1.066624e+04,  1.144149e-06},
    {1.00000e+06, 550.00, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP,  1.061570e+04,  1.110647e-06},
    {1.00000e+07, 550.00, {0.000000e+00, 0.000000e+00, 1.000000e+00, 0.000000e+00}, VAP,  1.017885e+04,  8.033311e-07},
    {1.01325e+05, 283.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -4.438247e+04, -5.187075e-07},
    {1.00000e+06, 283.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -4.426586e+04, -5.208575e-07},
    {1.00000e+07, 283.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -4.307920e+04, -5.386561e-07},
    {1.01325e+05, 298.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -4.060294e+04, -4.948567e-07},
    {1.00000e+06, 298.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -4.048960e+04, -4.975484e-07},
    {1.00000e+07, 298.15, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ, -3.933020e+04, -5.194386e-07},
    {1.01325e+05, 550.00, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, VAP,  6.288284e+04,  8.574199e-06},
    {1.00000e+06, 550.00, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, VAP,  6.015485e+04,  1.126202e-05},
    {1.00000e+07, 550.00, {0.000000e+00, 0.000000e+00, 0.000000e+00, 1.000000e+00}, LIQ,  3.752013e+04,  2.732726e-07},
    {1.01325e+05, 283.15, {3.286000e-01, 3.319000e-01, 3.313000e-01, 8.300000e-03}, VAP, -5.286490e+02,  5.366763e-06},
    {1.01325e+05, 283.15, {1.040000e-02, 3.000000e-04, 2.100000e-03, 9.871000e-01}, LIQ, -4.395667e+04, -5.167514e-07},
    {1.00000e+06, 283.15, {3.121000e-01, 3.468000e-01, 3.400000e-01, 1.000000e-03}, VAP, -6.521012e+02,  4.842468e-06},
    {1.00000e+06, 283.15, {9.190000e-02, 3.600000e-03, 2.110000e-02, 8.835000e-01}, LIQ, -4.041707e+04, -5.013302e-07},
    {1.00000e+07, 283.15, {1.850000e-01, 4.662000e-01, 3.480000e-01, 7.000000e-04}, VAP, -1.481792e+03,  2.276903e-06},
    {1.00000e+07, 283.15, {3.097000e-01, 5.140000e-02, 1.599000e-01, 4.790000e-01}, LIQ, -2.534339e+04, -3.983758e-07},
    {1.01325e+05, 298.15, {3.253000e-01, 3.276000e-01, 3.272000e-01, 1.990000e-02}, VAP, -1.910796e+01,  5.184416e-06},
    {1.01325e+05, 298.15, {8.000000e-03, 4.000000e-04, 1.900000e-03, 9.898000e-01}, LIQ, -4.029301e+04, -4.930559e-07},
    {1.00000e+06, 298.15, {3.171000e-01, 3.431000e-01, 3.374000e-01, 2.500000e-03}, VAP, -1.524379e+02,  4.467927e-06},
    {1.00000e+06, 298.15, {7.240000e-02, 3.800000e-03, 1.880000e-02, 9.050000e-01}, LIQ, -3.758874e+04, -4.793543e-07},
    {1.00000e+07, 298.15, {2.117000e-01, 4.398000e-01, 3.470000e-01, 1.500000e-03}, VAP, -1.021400e+03,  2.300945e-06},
    {1.00000e+07, 298.15, {2.901000e-01, 5.140000e-02, 1.485000e-01, 5.101000e-01}, LIQ, -2.399549e+04, -3.749308e-07},
    {1.01325e+05, 550.00, {2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01}, VAP,  2.296457e+04,  2.935798e-06},
    {1.00000e+06, 550.00, {2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01}, VAP,  2.268036e+04,  2.878424e-06},
    {1.00000e+07, 550.00, {2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01}, VAP,  2.014437e+04,  1.974289e-06},
  } )
);

INSTANTIATE_TEST_SUITE_P(
  CompositionalEnthalpyTest, SoaveRedlichKwong_6,
  ::testing::ValuesIn< typename SoaveRedlichKwong_6::DataType >( {
    {1.01325e+05, 283.15, {2.478590e-01, 1.239310e-01, 8.588000e-03, 2.478630e-01, 1.858620e-01, 1.858970e-01}, VAP, -7.928696e+02,  1.221540e-05},
    {1.01325e+05, 283.15, {1.400000e-05, 0.000000e+00, 9.998400e-01, 0.000000e+00, 1.450000e-04, 0.000000e+00}, LIQ, -4.821482e+04, -2.408303e-07},
    {1.00000e+06, 283.15, {2.498020e-01, 1.249180e-01, 1.023000e-03, 2.498360e-01, 1.870430e-01, 1.873770e-01}, VAP, -1.364451e+03,  1.298410e-05},
    {1.00000e+06, 283.15, {1.350000e-04, 0.000000e+00, 9.985240e-01, 1.000000e-06, 1.340000e-03, 0.000000e+00}, LIQ, -4.813521e+04, -2.411439e-07},
    {1.00000e+07, 283.15, {1.606730e-01, 3.112760e-01, 3.509000e-03, 4.560230e-01, 5.605000e-02, 1.247000e-02}, VAP, -2.718582e+03,  4.043744e-06},
    {1.00000e+07, 283.15, {2.139420e-01, 2.509800e-02, 2.696610e-01, 1.092340e-01, 1.833080e-01, 1.987580e-01}, LIQ, -2.059418e+04, -2.359039e-07},
    {1.01325e+05, 298.15, {2.441010e-01, 1.220520e-01, 2.361500e-02, 2.441050e-01, 1.830470e-01, 1.830790e-01}, VAP, -5.567760e+01,  1.118477e-05},
    {1.01325e+05, 298.15, {1.500000e-05, 0.000000e+00, 9.998430e-01, 0.000000e+00, 1.410000e-04, 0.000000e+00}, LIQ, -4.701695e+04, -2.356138e-07},
    {1.00000e+06, 298.15, {2.493660e-01, 1.247020e-01, 2.754000e-03, 2.494030e-01, 1.867230e-01, 1.870520e-01}, VAP, -5.727211e+02,  1.151918e-05},
    {1.00000e+06, 298.15, {1.500000e-04, 0.000000e+00, 9.985160e-01, 2.000000e-06, 1.333000e-03, 0.000000e+00}, LIQ, -4.693758e+04, -2.359289e-07},
    {1.00000e+07, 298.15, {1.885930e-01, 2.808340e-01, 6.328000e-03, 4.344900e-01, 7.191700e-02, 1.783800e-02}, VAP, -2.248771e+03,  4.125190e-06},
    {1.00000e+07, 298.15, {2.048090e-01, 2.376400e-02, 2.816490e-01, 1.011430e-01, 1.829180e-01, 2.057170e-01}, LIQ, -1.967401e+04, -1.819335e-07},
    {1.01325e+05, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  1.424903e+04,  3.202080e-06},
    {1.00000e+06, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  1.405756e+04,  3.173047e-06},
    {1.00000e+07, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  1.216120e+04,  2.686232e-06},
  } )
);

// This data is in mass units
INSTANTIATE_TEST_SUITE_P(
  CompositionalEnthalpyTest, SoreideWhitson_6,
  ::testing::ValuesIn< typename SoreideWhitson_6::DataType >( {  
    {1.01325e+05, 283.15, {2.478590e-01, 1.239310e-01, 8.588000e-03, 2.478630e-01, 1.858620e-01, 1.858970e-01}, VAP, -2.076518e+04,  1.253196e-05},
    {1.01325e+05, 283.15, {1.400000e-05, 0.000000e+00, 9.998400e-01, 0.000000e+00, 1.450000e-04, 0.000000e+00}, LIQ, -2.598884e+06, -2.346874e-07},
    {1.00000e+06, 283.15, {2.498020e-01, 1.249180e-01, 1.023000e-03, 2.498360e-01, 1.870430e-01, 1.873770e-01}, VAP, -3.590796e+04,  1.327629e-05},
    {1.00000e+06, 283.15, {1.350000e-04, 0.000000e+00, 9.985240e-01, 1.000000e-06, 1.340000e-03, 0.000000e+00}, LIQ, -2.591487e+06, -2.349740e-07},
    {1.00000e+07, 283.15, {1.606730e-01, 3.112760e-01, 3.509000e-03, 4.560230e-01, 5.605000e-02, 1.247000e-02}, VAP, -1.089422e+05,  4.162895e-06},
    {1.00000e+07, 283.15, {2.139420e-01, 2.509800e-02, 2.696610e-01, 1.092340e-01, 1.833080e-01, 1.987580e-01}, LIQ, -5.441170e+05, -2.191882e-07},
    {1.01325e+05, 298.15, {2.441010e-01, 1.220520e-01, 2.361500e-02, 2.441050e-01, 1.830470e-01, 1.830790e-01}, VAP, -1.510413e+03,  1.151375e-05},
    {1.01325e+05, 298.15, {1.500000e-05, 0.000000e+00, 9.998430e-01, 0.000000e+00, 1.410000e-04, 0.000000e+00}, LIQ, -2.537116e+06, -2.293975e-07},
    {1.00000e+06, 298.15, {2.493660e-01, 1.247020e-01, 2.754000e-03, 2.494030e-01, 1.867230e-01, 1.870520e-01}, VAP, -1.534557e+04,  1.183007e-05},
    {1.00000e+06, 298.15, {1.500000e-04, 0.000000e+00, 9.985160e-01, 2.000000e-06, 1.333000e-03, 0.000000e+00}, LIQ, -2.529758e+06, -2.296823e-07},
    {1.00000e+07, 298.15, {1.885930e-01, 2.808340e-01, 6.328000e-03, 4.344900e-01, 7.191700e-02, 1.783800e-02}, VAP, -8.762014e+04,  4.252175e-06},
    {1.00000e+07, 298.15, {2.048090e-01, 2.376400e-02, 2.816490e-01, 1.011430e-01, 1.829180e-01, 2.057170e-01}, LIQ, -5.183653e+05, -1.665260e-07},
    {1.01325e+05, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  4.147675e+05,  3.509469e-06},
    {1.00000e+06, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  4.086661e+05,  3.469178e-06},
    {1.00000e+07, 550.00, {2.000000e-01, 1.000000e-01, 2.000000e-01, 2.000000e-01, 1.500000e-01, 1.500000e-01}, VAP,  3.491882e+05,  2.852633e-06},
  } )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos

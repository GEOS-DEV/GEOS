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
  real64 const        // expected enthalpy
  >;

template< int NC >
struct FluidData {};

template<>
struct FluidData< 4 >
{
  static std::unique_ptr< TestFluid< 4 > > createFluid()
  {
    return TestFluid< 4 >::create( {Fluid::CO2, Fluid::H2, Fluid::CH4, Fluid::C2H6} );
  }

  static void populateCoefficients( HeatCapacityCoefficients * coefficients )
  {
    coefficients->m_referenceTemperature.resize( 4 );
    coefficients->m_referenceTemperature.zero();
    coefficients->m_referenceEnthalpy.resize( 1, 4 );
    coefficients->m_referenceEnthalpy.zero();
    coefficients->m_coefficients.resize( 1, 4, 5 );
    std::array< real64, 5*4 > coefficientsData{
      0.0, 0.0, 0.0, 0.0, 0.0,
      2.883, 0.003681, -7.720e-06, 6.920e-09, -2.130e-12,
      4.568, -0.008975, 3.631e-05, -3.407e-08, 1.091e-11,
      5.409, 0.1781, -0.00006938, 0.000000008713, 0.0
    };
//      4.178, -0.004427, 5.660e-05, -6.651e-08, 2.487e-11
//5.409, 0.1781, -0.00006938, 0.000000008713
    for( int ic = 0; ic < 4; ++ic )
    {
      for( int j = 0; j < 5; ++j )
      {
        coefficients->m_coefficients( 0, ic, j ) = coefficientsData[ic*5 + j];
      }
    }
  }
};

template< int NC, EquationOfStateType EOS_TYPE >
class CompositionalEnthalpyTestFixture :  public ::testing::TestWithParam< EnthalpyData< NC > >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numComps = NC;
  static constexpr int numDofs = NC + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  CompositionalEnthalpyTestFixture()
    : m_fluid( FluidData< NC >::createFluid() )
  {
    ComponentProperties const & componentProperties = this->m_fluid->getComponentProperties();
    m_parameters = CompositionalEnthalpy::createParameters( std::make_unique< ModelParameters >() );

    auto equationOfState = const_cast< EquationOfState * >(m_parameters->get< EquationOfState >());
    string const eosName = EnumStrings< EquationOfStateType >::toString( EOS_TYPE );
    equationOfState->m_equationsOfStateNames.emplace_back( eosName );

    auto heatCapacityCoefficients = const_cast< HeatCapacityCoefficients * >(m_parameters->get< HeatCapacityCoefficients >());
    FluidData< NC >::populateCoefficients( heatCapacityCoefficients );

    string const name = GEOS_FMT( "PhaseEnthalpy{}{}", eosName, NC );
    m_enthalpy = std::make_unique< CompositionalEnthalpy >( name, componentProperties, 0, *m_parameters );
  }

  ~CompositionalEnthalpyTestFixture() = default;

  void testEnthalpyValues( EnthalpyData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    //real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));
    //real64 const expectedEnthalpy = std::get< 3 >( data );

    auto componentProperties = m_fluid->createKernelWrapper();
    auto kernelWrapper = m_enthalpy->createKernelWrapper();

    real64 enthalpy = 0.0;
    stackArray1d< real64, numDofs > tempDerivs( numDofs );

    for( real64 t = 250.00; t <= 799.00; t += 50.0 )
    {
      kernelWrapper.compute( componentProperties,
                             pressure,
                             t,
                             phaseComposition.toSliceConst(),
                             enthalpy,
                             tempDerivs.toSlice(),
                             false );
      std::cout << std::fixed << std::setprecision( 0 ) << t << " "
                << std::fixed << std::setprecision( 5 ) << enthalpy << " "
                << std::fixed << std::setprecision( 5 ) << tempDerivs[1] << " "
                << std::endl;
    }
    //checkRelativeError( enthalpy, expectedEnthalpy, relTol, absTol );
  }

  void testEnthalpyDerivatives( EnthalpyData< NC > const & data )
  {
    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > phaseComposition;
    TestFluid< NC >::createArray( phaseComposition, std::get< 2 >( data ));

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
                           false );
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
                             false );
      return displacedEnthalpy;
    } );

    // -- Temperature derivative
    real64 const dT = 1.0e-6 * temperature;
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
                             false );
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
                               false );
        phaseComposition[ic] = z_old;
        return displacedEnthalpy;
      } );
    }
  }

protected:
  std::unique_ptr< CompositionalEnthalpy > m_enthalpy{};
  std::unique_ptr< TestFluid< NC > > m_fluid{};
  std::unique_ptr< ModelParameters > m_parameters{};
};

using PengRobinson = CompositionalEnthalpyTestFixture< 4, EquationOfStateType::PengRobinson >;

TEST_P( PengRobinson, testEnthalpyValues )
{
  testEnthalpyValues( GetParam() );
}
/**
   TEST_P( PengRobinson, testEnthalpyDerivatives )
   {
   testEnthalpyDerivatives( GetParam() );
   }
 **/
/* UNCRUSTIFY-OFF */

// Test data

INSTANTIATE_TEST_SUITE_P(
  CompositionalEnthalpyTest, PengRobinson,
  ::testing::ValuesIn<EnthalpyData<4>>( {
    {1.0e+05, 288.15, {0.000, 0.000, 0.000, 1.000}, 5.54544e+04}
  } )
);

/* UNCRUSTIFY-ON */

} // testing

} // geos

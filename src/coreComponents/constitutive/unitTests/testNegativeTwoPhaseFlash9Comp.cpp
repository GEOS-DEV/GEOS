/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 TotalEnergies
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

// Source includes
#include "codingUtilities/UnitTestUtilities.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 9;

using FlashData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  Feed< numComps > const,   // total composition
  bool,                     // expected flash status (success/failure)
  real64 const,             // expected vapour fraction
  Feed< numComps > const,   // expected liquid composition
  Feed< numComps > const    // expected vapour composition
  >;

template< typename EOS_TYPE >
class NegativeTwoPhaseFlashTest9CompFixture :  public ::testing::TestWithParam< FlashData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numDofs = numComps + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  NegativeTwoPhaseFlashTest9CompFixture()
    : m_fluid( createFluid() )
  {}

  ~NegativeTwoPhaseFlashTest9CompFixture() = default;

  void testFlash( FlashData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

    bool const expectedStatus = std::get< 3 >( data );
    real64 const expectedVapourFraction = std::get< 4 >( data );

    stackArray1d< real64, numComps > expectedLiquidComposition;
    stackArray1d< real64, numComps > expectedVapourComposition;
    TestFluid< numComps >::createArray( expectedLiquidComposition, std::get< 5 >( data ));
    TestFluid< numComps >::createArray( expectedVapourComposition, std::get< 6 >( data ));

    real64 vapourFraction = -1.0;
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray1d< real64, numComps > vapourComposition( numComps );
    stackArray2d< real64, numComps > kValues( 1, numComps );
    kValues.zero();

    bool status = NegativeTwoPhaseFlash::compute< EOS_TYPE, EOS_TYPE >(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition.toSlice(),
      vapourComposition.toSlice() );

    // Check the flash success result
    ASSERT_EQ( expectedStatus, status );

    std::cout << "EXPECTED " << vapourFraction << std::endl;

    if( !expectedStatus )
    {
      return;
    }
    GEOS_UNUSED_VAR( expectedVapourFraction );
    GEOS_UNUSED_VAR( expectedLiquidComposition );
    GEOS_UNUSED_VAR( expectedVapourComposition );
/**
    // Check the vaopur fraction
    checkRelativeError( expectedVapourFraction, vapourFraction, relTol, absTol );

    // Check liquid composition
    if( expectedVapourFraction < 1.0 - absTol )
    {
      for( integer ic=0; ic<numComps; ++ic )
      {
        checkRelativeError( expectedLiquidComposition[ic], liquidComposition[ic], relTol, absTol );
      }
    }

    // Check vapour composition
    if( absTol < expectedVapourFraction )
    {
      for( integer ic=0; ic<numComps; ++ic )
      {
        checkRelativeError( expectedVapourComposition[ic], vapourComposition[ic], relTol, absTol );
      }
    }
 */
  }

/**
   void testFlashDerivatives( FlashData< NC > const & data )
   {
    // Number of output values from each flash calculation
    constexpr integer numValues = 1 + 2*numComps;

    auto componentProperties = this->m_fluid->createKernelWrapper();

    bool const expectedStatus = std::get< 3 >( data );
    if( !expectedStatus ) return;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< NC >::createArray( composition, std::get< 2 >( data ));

    real64 vapourFraction = -1.0;
    stackArray1d< real64, numComps > liquidComposition( numComps );
    stackArray1d< real64, numComps > vapourComposition( numComps );

    stackArray1d< real64, numDofs > vapourFractionDerivs( numDofs );
    stackArray2d< real64, numComps * numDofs > liquidCompositionDerivs( numComps, numDofs );
    stackArray2d< real64, numComps * numDofs > vapourCompositionDerivs( numComps, numDofs );
    stackArray1d< real64, numValues > derivatives( numValues );

    // Combine values and derivatives into a single output
    auto const concatDerivatives = []( integer const kc, auto & derivs, auto const & v, auto const & xmf, auto const & ymf ){
      derivs[0] = v[kc];
      for( integer ic = 0; ic < numComps; ++ic )
      {
        derivs[1+ic] = xmf( ic, kc );
        derivs[1+ic+numComps] = ymf( ic, kc );
      }
    };
    std::cout << std::scientific << std::setprecision( 8 );

    auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ){
      stackArray1d< real64, numComps > displacedLiquidComposition( numComps );
      stackArray1d< real64, numComps > displacedVapourComposition( numComps );

      NegativeTwoPhaseFlash::compute< EOS_TYPE, EOS_TYPE >(
        numComps,
        p,
        t,
        zmf,
        componentProperties,
        values[0],
        displacedLiquidComposition,
        displacedVapourComposition );
      for( integer ic = 0; ic < numComps; ++ic )
      {
        values[1+ic] = displacedLiquidComposition[ic];
        values[1+ic+numComps] = displacedVapourComposition[ic];
      }
    };

    NegativeTwoPhaseFlash::compute< EOS_TYPE, EOS_TYPE >(
      numComps,
      pressure,
      temperature,
      composition,
      componentProperties,
      vapourFraction,
      liquidComposition,
      vapourComposition );

    NegativeTwoPhaseFlash::computeDerivatives< EOS_TYPE, EOS_TYPE >(
      numComps,
      pressure,
      temperature,
      composition,
      componentProperties,
      vapourFraction,
      liquidComposition,
      vapourComposition,
      vapourFractionDerivs,
      liquidCompositionDerivs,
      vapourCompositionDerivs );

    // Test against numerically calculated values
    // --- Pressure derivatives ---
    concatDerivatives( Deriv::dP, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
    real64 const dp = 1.0e-4 * pressure;
    geos::testing::internal::testNumericalDerivative< numValues >(
      pressure, dp, derivatives,
      [&]( real64 const p, auto & values ) {
      evaluateFlash( p, temperature, composition, values );
    } );

    // --- Temperature derivatives ---
    concatDerivatives( Deriv::dT, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
    real64 const dT = 1.0e-6 * temperature;
    geos::testing::internal::testNumericalDerivative< numValues >(
      temperature, dT, derivatives,
      [&]( real64 const t, auto & values ) {
      evaluateFlash( pressure, t, composition, values );
    } );

    // --- Composition derivatives ---
    real64 constexpr dz = 1.0e-7;
    for( integer jc = 0; jc < numComps; ++jc )
    {
      if( composition[jc] < 1.0e-6 ) continue;
      integer const kc = Deriv::dC + jc;
      concatDerivatives( kc, derivatives, vapourFractionDerivs, liquidCompositionDerivs, vapourCompositionDerivs );
      geos::testing::internal::testNumericalDerivative< numValues >(
        0.0, dz, derivatives,
        [&]( real64 const z, auto & values ) {
        real64 const originalFraction = composition[jc];
        composition[jc] += z;
        evaluateFlash( pressure, temperature, composition, values );
        composition[jc] = originalFraction;
      }, relTol, absTol );
    }
   }
 **/

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
private:
  static std::unique_ptr< TestFluid< numComps > > createFluid();
};

template< typename EOS_TYPE >
std::unique_ptr< TestFluid< numComps > > NegativeTwoPhaseFlashTest9CompFixture< EOS_TYPE >::createFluid()
{
  std::unique_ptr< TestFluid< numComps > > fluid = TestFluid< numComps >::create( {0, 0, 0, 0, 0, 0, 0, 0, 0} );
  // Manually populate
  TestFluid< numComps >::populateArray( fluid->criticalPressure, Feed< 9 >{73.8659e5, 33.9439e5, 46.0421e5, 48.8387e5, 42.4552e5, 37.47e5, 33.5892e5, 30.1037e5, 20.549e5} );
  TestFluid< numComps >::populateArray( fluid->criticalTemperature, Feed< 9 >{304.7, 126.2, 190.6, 305.43, 369.8, 419.5, 465.9, 507.5, 678.8} );
  TestFluid< numComps >::populateArray( fluid->criticalVolume, Feed< 9 >{9.3999e-05, 9.0001e-05, 9.7999e-05, 1.4800e-04, 2.0000e-04, 2.5800e-04, 3.1000e-04, 3.5100e-04, 6.8243e-04} );
  TestFluid< numComps >::populateArray( fluid->acentricFactor, Feed< 9 >{0.225, 0.04, 0.013, 0.0986, 0.1524, 0.1956, 0.2413, 0.299, 0.5618} );
  TestFluid< numComps >::populateArray( fluid->molecularWeight, Feed< 9 >{44.01e-3, 28.01e-3, 16.04e-3, 30.07e-3, 44.1e-3, 58.12e-3, 72.15e-3, 84e-3, 173e-3} );
  TestFluid< numComps >::populateArray( fluid->volumeShift, Feed< 9 >{ -0.04958, -0.136012, -0.1486264, -0.10863408, -0.08349872, -0.06331568, -0.04196464, -0.0150072, 0.0000 } );
  fluid->setBinaryCoefficients( Feed< 36 >{
        1.0000e-02,
        0.0000e+00, 3.7320e-03,
        0.0000e+00, 1.0000e-02, 0.0000e+00,
        0.0000e+00, 1.0000e-02, 0.0000e+00, 0.0000e+00,
        0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
        0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
        1.0000e-02, 0.0000e+00, 2.8000e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00,
        1.0000e-02, 0.0000e+00, 4.5320e-02, 1.0000e-02, 1.0000e-02, 0.0000e+00, 0.0000e+00, 0.0000e+00
      } );
  return std::move( fluid );
}

using NegativeTwoPhaseFlash9CompPR = NegativeTwoPhaseFlashTest9CompFixture< CubicEOSPhaseModel< PengRobinsonEOS > >;
using NegativeTwoPhaseFlash9CompSRK = NegativeTwoPhaseFlashTest9CompFixture< CubicEOSPhaseModel< SoaveRedlichKwongEOS > >;

TEST_P( NegativeTwoPhaseFlash9CompPR, testNegativeFlash )
{
  testFlash( GetParam() );
}

/**
   TEST_P( NegativeTwoPhaseFlash9CompSRK, testNegativeFlash )
   {
   testFlash( GetParam() );
   }
 */

std::vector< FlashData > generateTestData()
{
  std::vector< FlashData > data;
  auto pressures = { 1.0, 1.01325, 10.0, 100.0, 150.0, 500.0, 800.0, 1000.0 };
  auto temperatures = { 5.0, 15.5, 25.0, 50.0, 80.0, 100.0, 120.0, 250.0, 300.0, 600.0 };
  //auto pressures = { 340.0 };
  //auto temperatures = { 80.0 };
  const Feed< 9 > comp{ 0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920 };
  for( const auto & t : temperatures )
  {
    for( const auto & p : pressures )
    {
      data.emplace_back( FlashData{
            1.0e5*p, t + 273.15,
            comp,
            true, 0.0,
            comp, comp
          } );
    }
  }
  return data;
}


//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------
INSTANTIATE_TEST_SUITE_P(
  NegativeTwoPhaseFlash,
  NegativeTwoPhaseFlash9CompPR,
  ::testing::ValuesIn( generateTestData() )
  );

/**
   INSTANTIATE_TEST_SUITE_P(
   NegativeTwoPhaseFlash,
   NegativeTwoPhaseFlash9CompSRK,
   ::testing::Values(
    FlashData(
        1.000000e+05, 1.231500e+02,
        { 0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920 },
        true, 0.89038113,
        { 0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920 },
        { 0.00900, 0.00300, 0.53470, 0.11460, 0.08790, 0.04560, 0.02090, 0.01510, 0.16920 } )
    )
   );
 */

} // testing

} // geos

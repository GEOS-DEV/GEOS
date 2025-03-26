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
#include "constitutive/fluid/multifluid/compositional/functions/NegativeTwoPhaseFlash.hpp"
#include "constitutive/fluid/multifluid/compositional/functions/CubicEOSPhaseModel.hpp"
#include "TestFluid.hpp"
#include "TestFluidUtilities.hpp"

using namespace geos::constitutive::compositional;

namespace geos
{
namespace testing
{

static constexpr integer numComps = 3;

using FlashData = std::tuple<
  real64 const,             // pressure
  real64 const,             // temperature
  real64 const,             // salinity
  Feed< numComps > const,   // total composition
  real64 const,             // expected vapour fraction
  Feed< numComps > const,   // expected liquid composition
  Feed< numComps > const    // expected vapour composition
  >;

class SoreideWhitsonFlashTestFixture : public ::testing::TestWithParam< FlashData >
{
  static constexpr real64 relTol = 1.0e-5;
  static constexpr real64 absTol = 1.0e-7;
  static constexpr int numDofs = numComps + 2;
  using Deriv = geos::constitutive::multifluid::DerivativeOffset;
public:
  SoreideWhitsonFlashTestFixture()
    : m_fluid( TestFluid< numComps >::create( {Fluid::C1, Fluid::CO2, Fluid::H2O} ) )
  {
// These correlations only work with the correct values
    m_fluid->criticalPressure[0] = 46.0e5;
    m_fluid->criticalTemperature[0] = 190.6;
    m_fluid->acentricFactor[0] = 0.0108;
  }
  ~SoreideWhitsonFlashTestFixture() = default;

  void testFlash( FlashData const & data )
  {
    auto componentProperties = this->m_fluid->createKernelWrapper();

    constitutive::compositional::FlashData flashData;
    flashData.liquidEos = EquationOfStateType::SoreideWhitson;
    flashData.vapourEos = EquationOfStateType::PengRobinson;

    real64 const pressure = std::get< 0 >( data );
    real64 const temperature = std::get< 1 >( data );
    flashData.salinity = std::get< 2 >( data );
    stackArray1d< real64, numComps > composition;
    TestFluid< numComps >::createArray( composition, std::get< 3 >( data ));

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

    bool status = NegativeTwoPhaseFlash::compute(
      numComps,
      pressure,
      temperature,
      composition.toSliceConst(),
      componentProperties,
      flashData,
      kValues.toSlice(),
      vapourFraction,
      liquidComposition.toSlice(),
      vapourComposition.toSlice() );

    // Check the flash success result
    ASSERT_TRUE( status );

    if( !status )
    {
      return;
    }

    ((void)expectedVapourFraction);

    std::cout
      << "FlashData<" << numComps << ">{"
      << std::scientific << std::setprecision( 1 ) << pressure << ", "
      << std::fixed << std::setprecision( 2 ) << temperature << ", "
      << std::fixed << std::setprecision( 2 ) << flashData.salinity << ", "
      << std::fixed << std::setprecision( 3 ) << "{" << composition[0] << ", " << composition[1] << ", " << composition[2] << "}, "
      << std::fixed << std::setprecision( 8 ) << vapourFraction << ", "
      << std::fixed << std::setprecision( 6 ) << "{" << liquidComposition[0] << ", " << liquidComposition[1] << ", " << liquidComposition[2] << "}, "
      << std::fixed << std::setprecision( 6 ) << "{" << vapourComposition[0] << ", " << vapourComposition[1] << ", " << vapourComposition[2] << "}, "
                                                                                                                                              "},\n";

    // Check the vaopur fraction
    //checkRelativeError( expectedVapourFraction, vapourFraction, relTol, absTol );
/**
    // Check liquid composition
    if( expectedVapourFraction < 1.0 - absTol )
    {
      checkRelativeError( expectedLiquidComposition[0], liquidComposition[comp0], relTol, absTol );
      checkRelativeError( expectedLiquidComposition[1], liquidComposition[comp1], relTol, absTol );
    }

    // Check vapour composition
    if( absTol < expectedVapourFraction )
    {
      checkRelativeError( expectedVapourComposition[0], vapourComposition[comp0], relTol, absTol );
      checkRelativeError( expectedVapourComposition[1], vapourComposition[comp1], relTol, absTol );
    }
 **/
  }

  void testFlashDerivatives( FlashData const & data )
  {
    GEOS_UNUSED_VAR( data );
    /**
       // Number of output values from each flash calculation
       constexpr integer numValues = 1 + 2*numComps;

       auto componentProperties = this->m_fluid->createKernelWrapper();

       constitutive::compositional::FlashData flashData;
       flashData.liquidEos = EOS_TYPE;
       flashData.vapourEos = EOS_TYPE;

       bool const expectedStatus = std::get< 3 >( data );
       if( !expectedStatus ) return;

       real64 const pressure = std::get< 0 >( data );
       real64 const temperature = std::get< 1 >( data );
       stackArray1d< real64, numComps > composition;
       TestFluid< numComps >::createArray( composition, std::get< 2 >( data ));

       real64 vapourFraction = -1.0;
       stackArray1d< real64, numComps > liquidComposition( numComps );
       stackArray1d< real64, numComps > vapourComposition( numComps );
       stackArray2d< real64, numComps > kValues( 1, numComps );
       kValues.zero();

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

       auto const evaluateFlash = [&]( real64 const p, real64 const t, auto const & zmf, auto & values ){
       stackArray1d< real64, numComps > displacedLiquidComposition( numComps );
       stackArray1d< real64, numComps > displacedVapourComposition( numComps );

       NegativeTwoPhaseFlash::compute(
        numComps,
        p,
        t,
        zmf.toSliceConst(),
        componentProperties,
        flashData,
        kValues.toSlice(),
        values[0],
        displacedLiquidComposition.toSlice(),
        displacedVapourComposition.toSlice() );
       for( integer ic = 0; ic < numComps; ++ic )
       {
        values[1+ic] = displacedLiquidComposition[ic];
        values[1+ic+numComps] = displacedVapourComposition[ic];
       }
       };

       NegativeTwoPhaseFlash::compute(
       numComps,
       pressure,
       temperature,
       composition.toSliceConst(),
       componentProperties,
       flashData,
       kValues.toSlice(),
       vapourFraction,
       liquidComposition.toSlice(),
       vapourComposition.toSlice() );

       NegativeTwoPhaseFlash::computeDerivatives(
       numComps,
       pressure,
       temperature,
       composition.toSliceConst(),
       componentProperties,
       flashData,
       vapourFraction,
       liquidComposition.toSliceConst(),
       vapourComposition.toSliceConst(),
       vapourFractionDerivs.toSlice(),
       liquidCompositionDerivs.toSlice(),
       vapourCompositionDerivs.toSlice() );

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
       }, 10*relTol, 10*absTol );
       }*/
  }

protected:
  std::unique_ptr< TestFluid< numComps > > m_fluid{};
};

TEST_P( SoreideWhitsonFlashTestFixture, testNegativeFlash )
{
  testFlash( GetParam() );
}

std::vector< FlashData > generateTestData()
{
  std::vector< FlashData > testData;
  auto const feeds = {
    Feed< numComps >{0.0, 0.5, 0.5}
  };
  for( const auto & composition : feeds )
  {
    //for( const real64 pressure : {1.1e5, 1.0e6, 1.0e8} )
    real64 pressure = 1.0;
    real64 const dp = pow( 1000.0/pressure, 1.0/50 );
    for( int i = 0; i <= 50; i++ )
    {
      for( const real64 temperature : {2.97150e+02} )
      {
        for( const real64 salinity : {0.0} )
        {
          testData.emplace_back( 1e5*pressure, temperature, salinity, composition, 0.0, composition, composition );
        }
      }
      pressure *= dp;
    }
  }
  return testData;
}

//-------------------------------------------------------------------------------
// Data
//-------------------------------------------------------------------------------

/* UNCRUSTIFY-OFF */

INSTANTIATE_TEST_SUITE_P( SoreideWhitsonFlashTest, SoreideWhitsonFlashTestFixture, ::testing::ValuesIn( generateTestData()) );

/* UNCRUSTIFY-ON */

} // testing

} // geos
